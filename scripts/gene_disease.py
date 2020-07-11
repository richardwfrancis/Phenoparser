#!/usr/bin/env python3

import collections
#from collections import defaultdict
import requests
import requests_cache

requests_cache.uninstall_cache()
#from joblib import Parallel, delayed
#import multiprocessing
import time
#import sys
import json

gdata = collections.namedtuple('gdata', 'gene impact disease error')

class restClient(object):
    def __init__(self, reqs_per_sec=15):
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
        self.retry = 0
        self.retrymax = 3

    def perform_rest_action(self, url, params=None, resource='ensembl'):
        
        data = None
                
        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            #if delta < 1:
                #time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        
        r = None
        srch = ""
        if resource == "ensembl":
            r = requests.get(url, data=params, verify=True)
            srch = url
        else:
            r = requests.post(url, data=params, verify=True)
            srch = params['search']

        self.retry = self.retry + 1
            
        if r.ok:
            data = json.loads(r.text)
        else:
            # check if we are being rate limited by the server
            if r.status_code == 429:
                with open('acmg.log','a') as lf:
                    lf.write('\t\trate limited by {}\n'.format(resource))
                if 'Retry-After' in r.headers:
                    retry = r.headers['Retry-After']
                    with open('acmg.log','a') as lf:
                        lf.write('\t\tsleeping for {} for {}\n'.format(retry,resource))
                    time.sleep(float(retry))
                    with open('acmg.log','a') as lf:
                        lf.write('\t\tretrying {}\n'.format(resource))
                    data = self.perform_rest_action(url=url, params=params, resource=resource)                    
                else:
                    with open('acmg.log','a') as lf:
                        lf.write('\t\tno retry-after for {}\n'.format(srch))
                    data = json.loads( json.dumps({'reason':'no Retry-After for {} for {}'.format(resource,srch),'error':'{}'.format(r.status_code)}) )
            elif r.status_code == 400:
                dataerr = json.loads(r.text)
                if "error" in dataerr:
                    data = json.loads( json.dumps({'reason':dataerr["error"],'error':'{}'.format(r.status_code)}) )
                else:
                    if self.retry <= self.retrymax:
                        data = self.perform_rest_action(url=url, params=params, resource=resource)
                    else:
                        data = json.loads( json.dumps({'reason':'too many retries for {} for {}'.format(resource,srch),'error':'{}'.format(r.status_code)}) )
            else:
                data = json.loads( json.dumps({'reason':'unknown for {} for {}'.format(resource,srch),'error':'{}'.format(r.status_code)}) )

        return data


class geneClient(object):
    def __init__(self, gene=None):
        self.gene = gene
        self.gene_data_list = []
        self.impacts = []
        self.names = []
        

    def getSNPdata(self,impacts,snprs):
        ensembl_client = restClient()
        json_ensembl = ensembl_client.perform_rest_action(url='https://grch37.rest.ensembl.org/variation/human/'+snprs,
                                                          params={'content-type':'application/json'},
                                                          resource='ensembl')

        if "error" in json_ensembl:
            self.gene_data_list.append( gdata(gene=self.gene,impact="",disease="",error="ensembl failed with error code {} for {}".format(json_ensembl["error"], snprs)) )
        else:
            if 'most_severe_consequence' in json_ensembl:
                # if it's in there then process the contents
                self.impacts.append(json_ensembl['most_severe_consequence'])


    def parseAllelicVariant(self,value):
        for avar in value['allelicVariantList']:
            self.impacts = []
            self.names = []
            # find the interesting fields we want
            # make sure that this record is for the gene of interest
            # the only place I can find this is in the mutations string
            # it is usually the first item in a comma delimited list (e.g. "KCNQ1, ALA525THR")
            if 'mutations' in avar['allelicVariant']:
                mgene = avar['allelicVariant']['mutations'].split(',', 1)[0]
                if mgene == self.gene:
                    if 'name' in avar['allelicVariant']:
                        # if there are dbSNP IDS then head off to ensembl for the SNP data
                        if 'dbSnps' in avar['allelicVariant']:
                            # only append name if there are impacts to be found
                            self.names.append( avar['allelicVariant']['name'] )
                            # this is grch37
                            for snprs in avar['allelicVariant']['dbSnps'].split(','):
                                self.getSNPdata(self.impacts,snprs)
                else:
                    continue
            else:
                continue

            # associate SNP impact with disease for this variant
            for impact in self.impacts:
                for n in self.names:
                    self.gene_data_list.append( gdata(gene=self.gene,impact=impact,disease=n,error="") )              


    def getOmimData(self,omim_key):
        print("sleeping before starting work on {}".format(self.gene))
        time.sleep(2)
        print("waking up again to process {}".format(self.gene))
        if self.gene != 'None':
            omim_client = restClient()
            json_omim = omim_client.perform_rest_action(url='https://api.omim.org/api/entry/search?',
                                                        params={'search': self.gene,
                                                                'include': 'allelicVariantList',
                                                                'format':'json',
                                                                'apiKey': omim_key,
                                                                'limit': 5},
                                                        resource='omim')
            if "error" in json_omim:
                self.gene_data_list.append( gdata(gene=self.gene,impact="",disease="",error="omim failed with error code {} for {}".format( json_omim["error"], self.gene )) )
            else:
                for entry in json_omim['omim']['searchResponse']['entryList']:
                    for value in entry.values():
                        if 'allelicVariantList' in value:
                            self.parseAllelicVariant(value)
                        
                
def getGeneList(gene,omim_key):
    data = geneClient(gene=gene)
    # get the gene impact disease data from OMIM (only first 5 records returned for each gene - target gene should be the first)
    data.getOmimData(omim_key)
    print("completed {}".format(gene))
    return data.gene_data_list



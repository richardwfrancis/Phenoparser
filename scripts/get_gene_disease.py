#!/usr/bin/env python3

import sys
import getopt
from gene_disease import restClient, geneClient, getGeneList
import multiprocessing as mp
from nested_dict import nested_dict
import sqlite3
from sqlite3 import Error
import collections
import dill

def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
 
    return None

def make_table(conn):
    """
    Make the table to store the impact:disease data
    :param conn: the Connection object
    :return:
    """
    cur = conn.cursor()
    cur.execute('''create table if not exists impact_disease( gene TEXT default NULL, impact TEXT default NULL, disease text default NULL )''')
    cur.execute('''delete from impact_disease''')
 
    return None

def get_genes(conn,gdbtype):
    """
    Get the genes out of the gemini database or file
    :param conn: the Connection object or filename
    :param gdbtype: the type (gemini or file)
    :return:
    """
    allrows = []
    if gdbtype == "gemini":
        cur = conn.cursor()
        cur.execute("select distinct gene from gene_summary where is_hgnc = 1 and hgnc_id is not null")
        rows = cur.fetchall()
        for r in rows:
            allrows.append(r[0])
    else:
        f=open(conn,"r")
        allrows = [x.strip() for x in f.readlines()]
        f.close()

    return allrows

def load_table(conn,row_data):
    """
    Load the retrieved data
    :param conn: the Connection object
    :return:
    """
    cur = conn.cursor()

    for rowd in row_data:
        for row in rowd:
            if row.error == "":
                cur.execute('''insert into impact_disease(gene,impact,disease) VALUES (?,?,?)''', (row.gene,row.impact,row.disease))
            else:
                print("error from {} => {}".format(row.gene,row.error))

    conn.commit()

def getginfoasync(geneList,omim_key,cpu):
    g_information = []
    #pool = mp.Pool(mp.cpu_count())
    #pool = mp.Pool(100)
    pool = mp.Pool(int(cpu))
    #result_objects = [pool.apply_async(getGeneList, args=(gene, omim_key)) for gene in geneList]
    g_information = pool.starmap_async(getGeneList, [(gene, omim_key) for i, gene in enumerate(geneList)]).get()
    #print(results)
    # result_objects is a list of pool.ApplyResult objects
#    g_information = [r.get()[0] for r in result_objects]
    pool.close()
    #pool.join()

    return g_information

def getginfo(geneList,omim_key):
    g_information = []
    pool = mp.Pool(mp.cpu_count())
    #pool = mp.Pool(15)
    g_information = [pool.apply(getGeneList, args=(gene,omim_key)) for gene in geneList]
    #g_information = pool.starmap(getGeneList, [(gene, omim_key) for gene in geneList])
    #g_information = Parallel(n_jobs=100)(delayed(getGeneList)(gene,omim_key) for gene in geneList)
    pool.close()

    return g_information

def main(argv):
    gene_db = ""
    gene_db_type = ""
    impact_disease_db = "acmg_gid.db"
    omim_key = ""
    cpu = 2
    geneList = []

    try:
        opts, args = getopt.getopt(argv,"hk:g:t:o:c:",["omimkey=","genes=","gtype=","giddb=","cpu="])
    except getopt.GetoptError:
        print('get_gene_disease.py -k <omimKey> -g <geneData> -t <geneDataType> -o <geneImpactDiseaseDB>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('get_gene_disease.py -k <omimKey> -g <geneData> -t <geneDataType> -o <geneImpactDiseaseDB> -c <numThreads>')
            sys.exit()
        elif opt in ("-k", "--omimkey"):
            omim_key = arg
        elif opt in ("-g", "--genes"):
            gene_db = arg
        elif opt in ("-t", "--gtype"):
            gene_db_type = arg
        elif opt in ("-o", "--giddb"):
            impact_disease_db = arg
        elif opt in ("-c", "--cpu"):
            cpu = arg

    if ((omim_key == "" or gene_db == "" or gene_db_type == "")):
        print("You must provide both an omim key (-k) and a source of genes (-g).\nThis could be a GEMINI database as made by build_gemini.sh (-t gemini) or a file containing a list of genes (-t file)\nget_gene_disease.py -k <omimKey> -g <geneData> -t <geneDataType> -o <geneImpactDiseaseDB> -c <numThreads>")
        sys.exit()
    else:
        # create database connections
        # connect to the gene source if it is a database
        gene_conn = gene_db
        if (gene_db_type == "gemini"):
            print("connecting to {}".format(gene_db),flush=True)
            gene_conn = create_connection(gene_db)
            
        print("connecting to {}".format(impact_disease_db),flush=True)
        id_conn = create_connection(impact_disease_db)
        
        # make the table
        print("making impact_disease database",flush=True)
        make_table(id_conn)

        # get the genes from the database or file
        print("getting genes",flush=True)
        #geneList = ['NPPA','IRF6','MT-ND4','VWA1']
        #geneList = ['NPPA']
        geneList = get_genes(gene_conn,gene_db_type)

        # get the information
        print("getting gene information",flush=True)
        #g_information = getginfo(geneList,omim_key)
        g_information = getginfoasync(geneList,omim_key,cpu)
        print("collected information for {} genes".format(len(g_information)),flush=True)

        # load the data
        print("loading gene information",flush=True)
        load_table(id_conn,g_information)

        print("closing connections",flush=True)
        if (gene_db_type == "gemini"):
            gene_conn.close()
        id_conn.close()

if __name__ == "__main__":
   main(sys.argv[1:])

#!/usr/bin/env bash

pwd=$(pwd)

usage() {
    cat <<EOF
usage: $0  -i <inputdir>  -l <logfile>
EOF
    exit 1;
}

main() {

    LOG_STEPS=
    INDIR=
    NUM_THREADS=6
    
    while getopts d:i:t:l: opt
    do
        case ${opt} in
             
            i) INDIR=${OPTARG};;
            t) NUM_THREADS=${OPTARG};;
            l) LOG_STEPS=${OPTARG};;
       
            *) usage;;
        esac
    done

    if [ "${INDIR}" = "" ]; then usage; fi
    
    make_clean_copy_of_input;

    # run vt on input files... 
    for file in $(find tmp -name *.vcf -type f); do
        if ! [[ $file =~ ".d.vcf" ]]; then
            echo "$file running";            
            pipeline $file;
        fi;
    done   

    #
    #    Extract phenotype information - store in flatfiles... 
    #
    if file_exists sample_info.txt; then 
        rm sample_info.txt
    fi
    
    for file in $(find tmp -name *.d.n.vcf.gz -type f); do
        sanity_check_vcf_files $file;
    done
    exit;    
    echo "DONE!!! Hurrah ";

}

make_clean_copy_of_input() {
    
    local WORKINGDIR=$pwd/tmp/ 
    step "Create working directory";
    try make_working_directory $WORKINGDIR
    next

    # copy all relevant files into tmp directory.
    # !!!!!rename!!!!! samples if name in vcf is Unknown!!! 
    
    step "Make local copies of vcf and annotation files"
    next
    for file in $(find $INDIR -name *.vcf -type f); do
   
        local basefilename=$(basename $file)
        local SAMPLENAME=

        step "Working with:  $basefilename";
        if file_exists $WORKINGDIR/$basefilename; then
            echo "Warning: $WORKINGDIR/$basefilename exists.";
        else
       
            try cp $file $WORKINGDIR
            try chmod 700 $WORKINGDIR/$basefilename
            local SAMPLENAME=
            get_sampleID_from_vcf_file SAMPLENAME $WORKINGDIR/$basefilename;
       
            if [ "$SAMPLENAME" = "Unknown" ]; then
                local BASEDIR=$(dirname "$file") 
                local BASEDIR=$(echo $BASEDIR | rev | cut -d'/' -f2- | rev)
                local DIRNAME=$(basename $BASEDIR)
                SAMPLENAME="$DIRNAME";
                echo "WARNING: Unknown sample! Will rename to: $SAMPLENAME ! ";

                echo "$SAMPLENAME" > $WORKINGDIR/samplenames.tmp;
                try $BCFTOOLS/bcftools reheader  -s $WORKINGDIR/samplenames.tmp $WORKINGDIR/$basefilename  > $WORKINGDIR/tmp.vcf
                try cp $WORKINGDIR/tmp.vcf $WORKINGDIR/$basefilename
                try rm $WORKINGDIR/samplenames.tmp
                try rm $WORKINGDIR/tmp.vcf
           
            else
                echo "$SAMPLENAME";
            fi
       
            BASEDIR=$(dirname "$file") 
            BASEDIR=$(echo $BASEDIR | rev | cut -d'/' -f2- | rev)

            echo "$SAMPLENAME";
            #
            # Look for omim 
            #
            if [ -f $BASEDIR/pheno/omim.txt ]; then
                newname="$SAMPLENAME"
                newname+="_omim.txt"
           
                try cp $BASEDIR/pheno/omim.txt $WORKINGDIR/$newname
                try chmod 700 $WORKINGDIR/$newname
            fi
            #
            # Look for hpo file...
            #
            if [ -f $BASEDIR/pheno/hpo.txt ]; then
                newname="$SAMPLENAME"
                newname+="_hpo.txt"
           
                try cp $BASEDIR/pheno/hpo.txt $WORKINGDIR/$newname
                try chmod 700 $WORKINGDIR/$newname
            fi
        fi
        next
    done
}

sanity_check_vcf_files() {
    if ! [ "$1" ]
    then
        echo "Sanity Check needs an input vcf file";
        return 1;
    fi
    local basefilename=$(basename $1)

    BASEDIR=$(dirname "$1") 
    BASEDIR=$(echo $BASEDIR | rev | cut -d'/' -f2- | rev)

    echo "SANITY: $1";
    local SAMPLENAME=

    get_sampleID_from_vcf_file SAMPLENAME $1;
    
    DIRNAME=$(basename $BASEDIR)

    SOURCE=$( $BCFTOOLS/bcftools view -h $1 |  grep "^##source=" | sed 's/^##source=//' | perl -pe 's/[\n,\t, ]+/_/g' );

    REF=$( $BCFTOOLS/bcftools view -h $1 | grep "^##reference=" | sed 's/^##reference=//' | perl -pe 's/[\n,\t, ]+/_/g' );

    printf "%s\t%s\t%s\t%s\n"  "$1" "$SAMPLENAME" "$SOURCE" "$REF" >> sample_info.txt;
}

pipeline() {

    if ! [ "$1" ]
    then
        echo "Pipeline function needs an input vcf file";
        return 1;
    fi

    local INNAME=;
    local WORKINGDIR=$pwd/tmp;

    vt_pipeline $1 INNAME $WORKINGDIR;

# RF check to see if this can be replaced by index_vcf $WORKINGDIR $INNAME
# may be an issue with overwriting files (bgzip -f -c options are in the function but not below)

    step "Indexing..."
    if file_exists $WORKINGDIR/$INNAME; then
        try $HTSLIB/bgzip $WORKINGDIR/$INNAME
        try $HTSLIB/tabix -p vcf -f $WORKINGDIR/$INNAME.gz;
        try $GRABIX/grabix index $WORKINGDIR/$INNAME.gz;
    else
        echo "$WORKINGDIR/$INNAME.gz exists".
    fi

    next

    step "Docker cleanup"
    try cleanup_docker
    next
}

make_working_directory() {
    if ! [ "$1" ]
    then
        echo "mkdir function needs an input";
        return 1;
    fi
    echo "$1";
    mkdir -p $1;
    chmod 700 $1;
    return 0;
}

main "$@";

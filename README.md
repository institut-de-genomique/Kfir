# Kfir

Kfir is a program which compare a set of reads with a reference genome and yields the subset belonging to it.
Kfir takes in account complexity of kmer word. The complexity could be dust or shannon entropy.
Shannon entropy is a way

#Running Kfir
        Usage : kfir -g <genome file> -r <reads file> -o <output prefix> {Options}
        !! Note : Arguments with * are required, the other are optionnal !!
#Option
the script needs the argument with '*' the other are optionnal.

        *     -g  <str> : reference genome/chromosome/scaffold
        *     -r  <str> : first file of reads - or list of files (.fq or .fq.gz)
              -p  <str> : second file of reads if paired-end files - or list of files (.fq or .fq.gz)
              -k  <int> : word size for the dictionary, default is 19"
              -n  <int> : number of words in common to validate a read as belonging to the genome, default is 1
              -o  <str> : prefix of the output files

          LOW COMPLEXITY TEST : choose one of the following options, default no complexity test
              -d        : add a dust complexity mask on the read files
              -t<float> : if missing, no Shannon entropy comparison, if -t default threshold is 1.59
              -h        : this help

#Output
The script will output a fastq file with the selected reads
#Dependencies
The script is written in c++, it uses the library boost (http://www.boost.org/)

#License
Kfir is distributed open-source under CeCILL FREE SOFTWARE LICENSE. Check out http://www.cecill.info/ for more information about the contents of this license.

#Contact
kfir [a] genoscope [.] cns [.] fr

#Installation
Download the git repository

make

kfir -h

or run kfir : kfir -g <genome file> -r <reads file> -o <output prefix> {Options}


#Acknowledgments
This work was financially supported by the Genoscope, Institut de Genomique, CEA and Agence Nationale de la Recherche (ANR), and France Génomique (ANR-10-INBS-09-08).

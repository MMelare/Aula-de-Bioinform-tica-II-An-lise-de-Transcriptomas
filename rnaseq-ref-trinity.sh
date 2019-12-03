#!/bin/bash

# # Primeiramente deve-se passar um argumento na linha de comando após o nome do script, esse primeiro argumento refere-se
# ao diretório de entrada, ou seja, o diretório que contém as reads processadas

input=$1

#  Em seguida é feito uma validação para ver se esse argumento foi passado ou não. "if" sigifica "SE" e "!" significa negação
# então, se não foi passado o argumento, imprima na tela a mensagem entre aspas. Se foi passado mais não é um diretório, imprima na tela a mensagem entre aspas.



if [ ! ${input} ]
then   
        echo "Missing input directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input directory ${input}"
                exit
        fi
fi


# Também é preciso passar um segundo argumento, o output - diretório para armazenar o resultado do processo de montagem

output=$2

# De maneira análoga ao argumento anterior é realizado uma validação para o argumento ouput.

if [ ! ${output} ]
then   
        echo "Missing output directory"
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "Wrong output directory ${output}"
                exit
        fi
fi

# Número de processadores a serem utilizados 
num_threads="15"

# Espaço reservado para a montagem
mem_gb="16G"

# Declaração das variáveis referentes ao diretório de saída, ai diretório trinity_GG_input e ao diretório trinity_GG_assembled
basedir_out="${output}"

aligned_out="${basedir_out}/trinity_GG_input"

mkdir -p ${aligned_out}

trinity_out="${basedir_out}/trinity_GG_assembled"

# Criando diretórios para as saídas dos programas que serão utilizados a seguir


mkdir -p ${trinity_out}

if [ ! -e "${aligned_out}/All.sorted.bam" ]; then
	echo -e "Collecting alignments ..."

# Criação de outra variável para os arquivos de alinhamento bam gerados
	bamfiles=()

	bamfiles=( $( find ${input} -name 'Aligned.out.sorted.bam' ) )

# Fusão dos arquivos bam gerados em um único arquivo

	samtools merge -f ${aligned_out}/All.sorted.bam ${bamfiles[*]}

	samtools sort --threads ${num_threads} ${aligned_out}/All.sorted.bam > ${aligned_out}/All.csorted.bam
fi

if [ ! -d ${trinity_out}/Trinity.timing ]; then
	
	echo -e "Assembling step (Trinity) ..."

# Montagem utilizando o Trinity e o alinhamento contra o genoma de referência.

	Trinity --KMER_SIZE 23 \ #Tamanho dos kmers a serem gerados a partir das reads
		--output ${trinity_out} \  # Diretório de saída dos dados a serem gerados
		--seqType fq \ #Tipo de sequencia: fastq
		--max_memory ${mem_gb} \  # Máximo de memória destinada ao processo
		--CPU ${num_threads} \ # Número de processadores
		--min_contig_length 300 \  # Tamanho minimo dos contigs
		--genome_guided_bam ${aligned_out}/All.csorted.bam \  # Arquivo de alinhamento no formato bam
		--genome_guided_max_intron 10000 \    # Tamanho maximo dos introns
		 > ${trinity_out}/Trinity.log.out.txt \
		2> ${trinity_out}/Trinity.log.err.txt
		
fi

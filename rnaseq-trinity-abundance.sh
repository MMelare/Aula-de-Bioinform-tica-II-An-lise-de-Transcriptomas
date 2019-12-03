#!/bin/bash

# #Deve-se passar como primeiro argumento na linha de comando o diretório contendo os arquivos de entrada no formato .fastq processados e renomeados pelo awk

input=$1

# Se naõ for passado diretório input, print então na tela a mensagem dizendo que está faltando o argumento, e se foi 
# passado mais é o diretório errado ou não contem o arquivo necessário, então print na tela dizendo que o argumento 
# passado é errado

if [ ! ${input} ]
then   
        echo "Missing input (renamed for Trinity) directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input (renamed for Trinity) directory ${input}"
                exit
        fi
fi

# Como segundo argumento deve-se passar o diretório que foi armazenado o arquivo de montagem
trinity_output=$2

# De maneira análoga ao argumento anterior, é realizada a validação do argumento "trinity_output"

if [ ! ${trinity_output} ]
then   
        echo "Missing Trinity output directory"
        exit
else   
        if [ ! -d ${trinity_output} ]
        then   
                echo "Wrong Trinity output directory ${trinity_output}"
                exit
        fi
fi

# output - Como terceiro argumento deve-se passar o diretório para armazenar o resultado da avaliação de abundância
output=$3

# Também é realizada a validação do argumento "output"

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

# Número de processadores
num_threads="8"

# Criou-se uma variável para o diretório abundance, criado dentro do diretório de saída

abundance_out="${output}/abundance"

# Criou-se o diretório abundance dentro do diretório de saída
mkdir -p ${abundance_out}

# Criou-se 2 variáveis "left" e "right" para o conjunto de reads1 e reads 2 respectivamente
left=()
right=()

echo "Collecting reads step ..."


left=($(find ${input} -type f -name '*.prinseq_1.fastq'))

# Remoção dos arquivos que não serão utilizados

rm -f ${abundance_out}/samples.txt
rm -f ${abundance_out}/quant_files.txt
rm -f ${abundance_out}/groups.txt

echo -e "id\tname\tgroup" > ${abundance_out}/groups.txt

# Para pegar o conjunto de reads 2:
for l in ${left[@]}; do
	repname=`basename ${l} | sed 's/\..*$//'`
	condname=`echo ${repname} | sed 's/[0-9]\+//'`
	r=`echo ${l} | sed 's/_1.fastq/_2.fastq/'`
	right=(${right[@]} ${r})

	echo -e "${condname}\t${abundance_out}/${repname}\t${l}\t${r}" >> ${abundance_out}/samples.txt
	echo -e "${abundance_out}/${repname}/quant.sf" >> ${abundance_out}/quant_files.txt

	echo -e "${repname}\t${repname}\t${condname}" >> ${abundance_out}/groups.txt
done


#echo ${left[*]}
#echo ${right[*]}

# Declaração de mais duas variáveis geradas pelo Trinity:

trinity_fasta=`find ${trinity_output} -type f -name 'Trinity*.fasta'`
trinity_trans_map=`find ${trinity_output} -type f -name Trinity*.gene_trans_map`

echo "Estimating abundances ..."

# Script para estimar a abundancia dos genes:

${TRINITY_HOME}/util/align_and_estimate_abundance.pl 	--transcripts	${trinity_fasta} \ # Arquivo de montagem
							--est_method	salmon \  # Método de estimação da abundancia
							--salmon_add_opts "--validateMappings" \ # Opções adicionais para o método de estimação escolhido
							--samples_file	${abundance_out}/samples.txt \ # Arquivo contendo as reads, explicitando as condições biológicas e réplicas.
							--gene_trans_map ${trinity_trans_map} \ #Arquivo gene_trans_map gerado na montagem pelo Trinity
							--prep_reference \  # Geração de um index do genoma de referencia
							--thread_count ${num_threads} \ # Número de processadores
							--seqType fq \ # Tipo de Sequencia: fastq
							--output_dir ${abundance_out} \ #Diretório de saída
							 > ${abundance_out}/align_and_estimate_abundance.log.out.txt \
							2> ${abundance_out}/align_and_estimate_abundance.log.err.txt

echo "Constructing abundance matrix ..."

# Para geração de uma matriz da abundancia estimada:

${TRINITY_HOME}/util/abundance_estimates_to_matrix.pl	--est_method salmon \ # Método de estimação
							--gene_trans_map ${trinity_trans_map} \ #Arquivo gene_trans_map gerado na montagem com o Trinity
							--name_sample_by_basedir \ As amostras são nomeadas de acordo com o diretório e não com o arquivo
							--cross_sample_norm none \ # Método a ser utilizado para a normalização dos genes diferencialmente expressos
							--quant_files ${abundance_out}/quant_files.txt \ # Diretório que irá receber os arquivos target gerados
							--out_prefix ${abundance_out}/abundance \ # Diretório que contem o arquivo de abundancia gerado para se fazer a matriz
							 > ${abundance_out}/abundance_estimates_to_matrix.log.out.txt \
							2> ${abundance_out}/abundance_estimates_to_matrix.log.err.txt 
		
echo "Calculating Differentially Expressed Genes ..."

# Criação do diretório DEG dentro do diretório abundance
mkdir -p ${abundance_out}/DEG

run-DESeq2.R 	--in="${abundance_out}/abundance.gene.counts.matrix"  \ # Input - matriz de estimação da abundancia gerada
		--groups="${abundance_out}/groups.txt" \  #Agrupamentos gerados
		--out="${abundance_out}/DEG" \  #output
		 > ${abundance_out}/DEG/run-DESeq2.log.out.txt \
		2> ${abundance_out}/DEG/run-DESeq2.log.err.txt



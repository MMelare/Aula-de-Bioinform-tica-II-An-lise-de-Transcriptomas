#!/bin/bash

# # Primeiramente deve-se passar um argumento na linha de comando após o nome do script, esse primeiro argumento refere-se
# ao diretório de entrada, ou seja, o diretório que contém as reads processadas
input=$1

#  Em seguida é feito uma validação para ver se esse argumento foi passado ou não. "if" sigifica "SE" e "!" significa negação
# então, se não foi passado o argumento, imprima na tela a mensagem entre aspas. Se foi passado mais não é um diretório,
# imprima na tela a mensagem entre aspas.

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

# De maneira análoga a validação do input, tem-se a validação do argumento "output"
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

# Criou-se uma variável para se referir ao número de processadore a ser utilizado
num_threads="8"

# Criou-se outra varável para se referir ao espaço destinado a montagem, quanto de memória pode ser destinado a esse processo.
mem_gb="10G"

# Criou-se outra variável para se referir ao diretório passado como output e outras 2 para 2 diretórios a serem criados dentro do ouput.

basedir_out="${output}/"

renamed_out="${basedir_out}/renamed"

trinity_out="${basedir_out}/trinity_assembled"

# Para criar o diretório renamed:
mkdir -p ${renamed_out}

# Criou-se outras 4 variáveis para se referir ao conjunto de reads pareadas e não pareadas 1 e 2.

left=()
left_singleton=()

right=()
right_singleton=()

echo "Renaming step ..."

# Criação do diretório trinity_assembled:

mkdir -p ${trinity_out}

# Etapa de renomeação das reads processadas. É preciso que essas estejam com o sufixo /1 e /2 para as readas 1 e 2 respectivamente.
# Para isso utilizou-se o comando abaixo "awk", em uma estrutura de loop, ou seja, para todos os arquivos que terminam com .fastq dentro do diretório de entrada, faça a renomeação.

for fastq in `ls ${input}/*.fastq`; do
	# obtendo nome do arquivo 
	fastqbn=`basename ${fastq}`;
	if [[ ! $fastqbn =~ \.bad_ ]]; then
		renamed_fastq="${renamed_out}/${fastqbn}"
		if [ ! -e ${renamed_fastq} ]; then
			echo -e "\tRenaming ${fastqbn} ..."
			if [[ ${fastqbn} =~ _1[\._] ]]; then
				awk '{ if (NR%4==1) { if ($1!~/\/1$/) { print $1"/1" } else { print $0 } } else if (NR%4==3) { print "+" } else { print $0 } }' ${fastq} > ${renamed_fastq}
			elif [[ ${fastqbn} =~ _2[\._]  ]]; then
				awk '{ if (NR%4==1) { if ($1!~/\/2$/) { print $1"/2" } else { print $0 } } else if (NR%4==3) { print "+" } else { print $0 } }' ${fastq} > ${renamed_fastq}
			else 
				echo "Warning: ${fastqbn} discarded!"
			fi
		fi
		
		if [[ ${fastqbn} =~ _1[\._] ]]; then
			if [[ ${fastqbn} =~ singletons ]]; then
				left_singleton=($(printf "%s\n" ${left_singleton[@]} ${renamed_fastq} | sort -u ))
			else
				left=($(printf "%s\n" ${left[@]} ${renamed_fastq}  | sort -u ))
			fi
		elif [[ ${fastqbn} =~ _2[\._] ]]; then
			if [[ ${fastqbn} =~ singleton ]]; then
				right_singleton=($(printf "%s\n" ${right_singleton[@]} ${renamed_fastq}  | sort -u ))
			else
				right=($(printf "%s\n" ${right[@]} ${renamed_fastq}  | sort -u ))
			fi
		else
			echo "Warning: ${fastqbn} discarded!"
		fi
	fi
done

# Se ainda não existe o arquivo Trinity.timing realize o comando abaixo:

if [ ! -d ${trinity_out}/Trinity.timing ]; then
	
	echo -e "Assembling step (Trinity) ..."
	
	rm -fr ${trinity_out}/
	mkdir -p ${trinity_out}

# Montagem de novo utilizando o programa Trinity:
	
	Trinity --KMER_SIZE 27 \ # Tamanho do kmers a serem gerados a partir das reads
		--output ${trinity_out} \  # Diretório de saída
		--seqType fq \  # Formato da sequencia: fastq
		--max_memory ${mem_gb} \  # Máximo de espaço disponivel
		--CPU ${num_threads} \  # Número de processadores
		--min_per_id_same_path 95 \ # Minimo de identidade para 2 sequencias serem unidas em uma única sequencia
		--max_diffs_same_path  5 \  # Maximo de diferenças permitidas para 2 sequencias se unirem
		--path_reinforcement_distance 5 \ # Minimo de sobreposição 
		--group_pairs_distance 500 \ # Distancia  para que as reads sejam tratadas como paired-end
		--min_glue 5 \ # Numero minimo de reads para que se unam os contigs
		--min_contig_length 600 \ #Tamanho minimo do contig
		--min_kmer_cov 3 \  # minimo de kmers para serem montados
		--left $(IFS=, ; echo "${left[*]},${left_singleton[*]}") \ # Reads 1
		--right $(IFS=, ; echo "${right[*]},${right_singleton[*]}") \ # Reads 2
		 > ${trinity_out}/Trinity.log.out.txt \
		2> ${trinity_out}/Trinity.log.err.txt
fi

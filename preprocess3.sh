#!/bin/bash

# Primeiramente deve-se passar um argumento na linha de comando após o nome do script, esse primeiro argumento refere-se
# ao diretório de entrada, ou seja, o diretório que contém as reads provenientes do sequenciamento. As mesmas serão processadas.

indir=$1

# Em seguida é feito uma validação para ver se esse argumento foi passado ou não. "if" sigifica "SE" e "!" significa negação
# então, se não foi passado o argumento, imprima na tela a mensagem entre aspas. Se foi passado mais não é um diretório,
# imprima na tela a mensagem entre aspas.

if [ ! ${indir} ]; then
	echo "Missing input directory."
	exit
fi


if [ ! -d ${indir} ]; then
	echo "Wrong input directory (${indir})."
	exit
fi

# Deve-se passar como segundo argumento na linha de comando de execução do script um diretório de saída, ou seja,
# onde serão armazenados os dados gerados.

outdir=$2

# De maneira análoga a validação anterior, o mesmo é feito para esse segundo argumento:

if [ ! ${outdir} ]; then
	echo "Missing output directory."
	exit
fi

# SE ${outdir} NÃO É DIRETÓRIO
if [ ! -d ${outdir} ]; then
	echo "Wrong output directory (${outdir})."
	exit
fi

# A função mkdir é responsável pela criação de diretórios. 
# Esses diretórios serão criados dentro do diretório passado como segundo argumento, como output, na linha de execução do script.
# Dentro do diretório passado, será criado outro chamado "processed" e dentro desse, serão criados 3 diretórios:
# FASTQC, ATROPOS e PRINSEQ. Sendo que no diretório fastqc é criado mais 2 diretórios: "pre" e "pos" referentes ao 
# pre e pos processamento.

mkdir -p ${outdir}/processed/fastqc/pre
mkdir -p ${outdir}/processed/atropos
mkdir -p ${outdir}/processed/prinseq
mkdir -p ${outdir}/processed/fastqc/pos

# O "for" é utilizado como um loop, ou seja, uma repetição, o seguinte comando diz para localizar dentro do diretório de
# entrada todos os arquivos que terminam com "_R1.fastq", e "_R2.fastq". Enquanto tiver esses arquivos continue procurando
# e fazendo a seguinte função.

for r1 in `ls ${indir}/*_R1.fastq`; do
	
	r2=`echo ${r1} | sed 's/_R1.fastq/_R2.fastq/'`
	
	if [ ! -e "${r2}" ]; then
		echo "Read 2 (${r2}) paired with Read 1 ($r1) not found."
		exit
	fi

	name=`basename ${r1} | sed 's/_R1.fastq//'`

# O comando echo serve para imprimir na tela qualquer mensagem, essa deve vim entre aspas.
	echo -e "FastQC pre-evaluation using sample ${name}: ${r1} & ${r2} ...\n"

# Aqui será realizada a etapa de checagem de qualidade utilizando o programa fastqc, essa primeira checagem é realizada antes
# da etapa de processamento.É feito tanto para as reads 1, quanto para as reads2.

	fastqc -t 2 \
   		${r1} \
   		-o ${outdir}/processed/fastqc/pre/ \
		 > ${outdir}/processed/fastqc/pre/${name}_R1.log.out.txt \
		2> ${outdir}/processed/fastqc/pre/${name}_R1.log.err.txt

	fastqc -t 2 \
   		${r2} \
		-o ${outdir}/processed/fastqc/pre/ \
		 > ${outdir}/processed/fastqc/pre/${name}_R2.log.out.txt \
		2> ${outdir}/processed/fastqc/pre/${name}_R2.log.err.txt
	
	echo -e "Running atropos (insert) for adapter trimming using sample ${name}: ${r1} & ${r2} ...\n"

# Em seguida, será realizado a trimagem dos adaptadores utilizando o software Atropos.
# Para ativá-lo é preciso ativar o ambiente pyenv, pois o programa foi desenvolvido em python;

	eval "$(pyenv init -)"

	pyenv activate atropos

# A primeira etapa do Atropos é a remoção dos adaptadores por meio do algoritmo do tipo 'insert match', onde
# as readas são alinhadas, ou seja, a read 1 é alinhadas com o reverso da read 2 e então o que for remanescente será trimado.

	atropos trim --aligner insert \
             -e 0.1 \  # Taxa de erro maxima permitida durante o alinhamento
             -n 2 \   # Número máximo de adaptadores para ser trimados por read
             -m 10 \   # Tamanho minimo das reads
             --op-order GAWCQ \ #Ordem de execução do processo
             --match-read-wildcards \  #Habilita interpretação da IUPAC: Aceita bases coringas
             -O 20 \  # Sobreposição minima das reads para considerar um alinhamento válido
             -q 27 \  # Score de qualidade minimo, abaixo disso bases das extremidades são trimadas
             -T 8 \ # Número de processadores
             --correct-mismatches conservative \ # Metodologia a ser utilizada com os mismatches encontrados no alinhamento
             --pair-filter any \ # Qual read será avaliada no processo de filtragem
             -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \ # Sequência de adaptadores da read1 na extremidade 3'
             -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT  \  # Sequencia de adaptadores da read 2 na extremidade 3'
             -o ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \ # output das reads1 que foram trimadas
             -p ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \  # output das reads2 que foram trimadas
             -pe1 ${r1} \  # Reads 1 
             -pe2 ${r2} \  # Reads 2 
             --untrimmed-output        ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq \ # Reads1 que não foram trimadas, servirão de entrada na próxima etapa;
             --untrimmed-paired-output ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq \ # Reads2 que	não foram trimadas, servirão de	entrada	na próxima etapa;
               > ${outdir}/processed/atropos/${name}.atropos.log.out.txt \
              2> ${outdir}/processed/atropos/${name}.atropos.log.err.txt
	
	echo -e "Running atropos (adapter) for adapter trimming using sample ${name}: ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq & ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq ...\n"

# Na segunda etapa, as reads que não forem trimadas são submetidas a trimagem do tipo 'adapter-match', onde, é realizado
# o alinhamento entre os adaptadores, por correspondencia, são entçao trimados.
	atropos trim    --aligner adapter \
                -e 0.1 \  # valor máximo de taxa de erro permitida no alinhamento
                -n 2 \    # Máximo de adaptadores removidos por read
                -m 10 \   # Tamanho minimo das reads para continuarem no processo. Abixo desse valor são excluidas.
                --match-read-wildcards \ # Modo de interpretação da IUPAC é ativado.
                -O 3 \    # Tamanho minimo de sobreposição das reads para que essas sejam consideradas alinhadas e em seguida trimadas.
                -q 20 \   # Valor de score de qualidade minimo. Abaixo desse valor as extremidades serão trimadas.
                --pair-filter both \  # Escolhe qual das reads entre na avaliação para serem filtradas.
                -pe1 ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq \ # Reads 1 que não foram trimadas
                -pe2 ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq \ # Reads 2 que não foram trimadas.
                -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \ # sequência do adaptador a ser removido na primeira read na região 3’
                -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT  \ # sequência do adaptador a ser removido na primeira read na região 3’
                -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \ # Sequência do adaptador a ser removido na read 2 na extremidade 3'.
                -G CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \ # Sequência do adaptador a ser removido na read 2 na extremidade 5'.
                -T 2 \ # Número de processadores 
                -o  ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq  \ # output das reads 1
                -p  ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq \  # output das reads 2
                 >  ${outdir}/processed/atropos/${name}.atropos_adapter.log.out.txt \
                2>  ${outdir}/processed/atropos/${name}.atropos_adapter.log.err.txt
	
	echo -e "Merging atropos adapter trimming results using sample ${name}: ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq and ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq + ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq and ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq ...\n"

# Fusão das reads 1 que foram trimadas em ambos os processos em um único arquivo.

	cat       ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \
        	  ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq \
   		> ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq

# Fusão	das reads 2 que foram trimadas em ambos os processos em um único arquivo.

	cat       ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \
        	  ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq \
   		> ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq

	echo -e "Removing useless atropos results ...\n"

# Remoção de arquivos que não serão utilizados

	rm -f ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \
	      ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq \
	      ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \
	      ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq

	echo -e "PrinSeq processing: ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq & ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq ...\n"

# Remoção de bases de baixa qualidade, adaptadores, artefatos utilizando o programa PrinSeq

	prinseq-lite.pl -fastq  ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq \ # Arquivos fastq trimados pelo Atropos - reads1
			-fastq2 ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq \ #Arquivos fastq trimados pelo atropos - reads 2
			-out_format 3 \ # Formato do output - fastq
			-trim_qual_window 3 \ # Janela a ser considerada para calculo do score de qualidade
			-trim_qual_step 1 \  # Qual o passo para calculo de qualidade (de quantas em quantas bases"
			-trim_qual_right 30 \ # Trima as bases com score abaixo desse valor na extremidade 3'.
			-trim_qual_type mean \ # Base para cálculo do score de qualidade: Média
			-trim_qual_rule lt \  #  Regra a ser seguida para comparação dos valores de score de qualidade
			-out_good ${outdir}/processed/prinseq/${name}.atropos_final.prinseq \ #Diretório de saída dos resultados positivos e nome do arquivo
			-out_bad  null \ # Diretório de saída dos resultados negativos
			-lc_method dust \ # Método a ser utilizado para filtragem de sequências de baixa complexidade.DUST:Os scores são computados com base em diferentes trinucleotideos.
			-lc_threshold 30 \ # Valor máximo para filtragem de sequencias
			-min_len 20 \ # Trimagem de sequencias menores que esse valor
			-trim_tail_right 5 \ # Bases a serem trimadas na extremidade 3'
			-trim_tail_left 5 \ # Bases a serem trimadas na extremidade 5'
			-noniupac \ # Exclui bases que não sejam A,T,C,G e N
			 > ${outdir}/processed/prinseq/${name}.atropos_final.prinseq.out.log \
			2> ${outdir}/processed/prinseq/${name}.atropos_final.prinseq.err.log

	echo -e "FastQC pos-evaluation using sample ${name}: ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq & ${outdir}/processed/prinseq/${name}_2.atropos_final.prinseq.fastq ...\n"

# Em seguida, é realizada novamente a etapa de checagem de qualidade das reads obtidas após o processamento. 
# As mesmas podem ser comparadas antes do processamento.

	fastqc -t 2 \
	   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq \
	   -o ${outdir}/processed/fastqc/pos/ \
	    > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1.log.out.txt \
	   2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1.log.err.txt

	fastqc -t 2 \
	   ${outdir}/processed/prinseq/${name}_2.atropos_final.prinseq.fastq \
	   -o ${outdir}/processed/fastqc/pos/ \
	    > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2.log.out.txt \
	   2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2.log.err.txt

# Se tiver reads que não possuem seu respectivo par - singletons - renomear

	# SE EXISTIR <SAMPLE_NAME>.atropos_final.prinseq_1_singletons.fastq
	if [ -e "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq" ]; then
		fastqc -t 2 \
		   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq \
	   	   -o ${outdir}/processed/fastqc/pos/ \
	            > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons.log.out.txt \
	           2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons.log.err.txt
	fi

	# SE EXISTIR <SAMPLE_NAME>.atropos_final.prinseq_2_singletons.fastq
	if [ -e "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq" ]; then
		fastqc -t 2 \
		   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq \
		   -o ${outdir}/processed/fastqc/pos/ \
	            > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2_singletons.log.out.txt \
	           2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2_singletons.log.err.txt
	fi


done

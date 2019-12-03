#!/bin/bash

# Primeira variável criada: Número de processadores a serem utilizados durante a execução do script.
# Quando se cria uma variável, facilita, e em todos os comandos que necessitar usá-la, deve-se colocar a variável em vez do número, ou string.

num_threads=4

# Na linha de comando de execução do script, é preciso passar um argumento que indica o diretório input, ou seja, as reads do sequenciamento

indir=$1


# Validação: SE(if) NÃO (!) FOI PASSADO ARGUMENTO 1 NA LINHA DE COMANDO, então imprima na tela a mensagem entre aspas.
if [ ! ${indir} ]; then
	echo "Missing input directory."
	exit
fi

# SE o argumento passado NÃO É um  DIRETÓRIO, então imprima na tela a mensagem entre aspas.
if [ ! -d ${indir} ]; then
	echo "Wrong input directory (${indir})."
	exit
fi

# A execução do script também requer um segundo argumento, o diretório de saída dos resultados gerados:

outdir=$2

# Da mesma maneira, é realizado uma validação para esse argumento:
# SE NÃO FOI PASSADO ARGUMENTO 2 NA LINHA DE COMANDO, entre imprima na tela a mensagem entre aspas.
if [ ! ${outdir} ]; then
	echo "Missing output directory."
	exit
fi

# SE o argumento passado não é um DIRETÓRIO, imprima na tela a mensagem entre aspas.
if [ ! -d ${outdir} ]; then
	echo "Wrong output directory (${outdir})."
	exit
fi

# Como terceiro argumento na linha de comando de execução do script, deve-se passar a anotação do genoma de referência no formato gtf:

refgtf=$3

# Novamente é realizado a validação para o terceiro argumento:
# SE NÃO FOI PASSADO ARGUMENTO 3 NA LINHA DE COMANDO, então imprima na tela a mensagem entre aspas.
if [ ! ${refgtf} ]; then
	echo "Missing GTF file."
	exit
fi

# Se o argumento não é um arquivo, então imprima na tela a mensagem entre aspas.
if [ ! -e "${refgtf}" ]; then
	echo "Not found GTF file (${refgtf})."
	exit
fi

# A linha de execução do script requer também um quarto argumento, o arquivo FASTA do Genoma de referência:

refseq=$4

# Validação:SE NÃO FOI PASSADO ARGUMENTO 4 NA LINHA DE COMANDO, imprima na tela a mensagem entre aspas.
if [ ! ${refseq} ]; then
	echo "Missing GENOME fasta file."
	exit
fi

# Se o arquivo não existe, então imprima na tela a mensagem entre aspas.
if [ ! -e "${refseq}" ]; then
	echo "Not found GENOME fasta file (${refseq})."
	exit
fi

# Esse script realiza uma chamada de outro script, no caso o script preprocess3.sh que é referente ao processamento dos dados.

./preprocess3.sh "${indir}" "${outdir}"

# A função mkdir cria diretórios
# Dentro do diretótio passado como output, como segundo argumento, serão criados os seguintes diretórios:

mkdir -p ${outdir}/star_index
mkdir -p ${outdir}/star_out_pe
mkdir -p ${outdir}/star_out_se
mkdir -p ${outdir}/star_out_final
mkdir -p ${outdir}/cufflinks
mkdir -p ${outdir}/cuffmerge
mkdir -p ${outdir}/stringtie
mkdir -p ${outdir}/stringmerge
mkdir -p ${outdir}/cuffcompare
mkdir -p ${outdir}/cuffquant

# O "for" é utilizado como uma estrtura de loop, facilitando para que nao se tenha q fazer manualmente toda a linha de comando, colocando todas as reads.
# Para a variável r1 busque no diretório prinseq todos os arquivos que finalizam com ".atropos_final.prinseq_1.fastq", e faça o comando seguinte.
# Para as reads singletons, substitua "prinseq_1.fastq" por "prinseq_1_singletons.fastq/".
# Se não existe singletons, então crie um arquivo vazio referente aos singletons.
# Para r2, faça o mesmo mais para as reads 2.

for r1 in `ls ${outdir}/processed/prinseq/*.atropos_final.prinseq_1.fastq`; do
	r1_singletons=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_1_singletons.fastq/'`
	if [ ! -e "${r1_singletons}" ]; then
		touch ${r1_singletons}
	fi
	
	r2=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_2.fastq/'`
	
	if [ ! -e "${r2}" ]; then
		echo "Read 2 (${r2}) paired with Read 1 ($r1) not found."
		exit
	fi
	
	r2_singletons=`echo ${r2} | sed 's/prinseq_2.fastq/prinseq_2_singletons.fastq/'`
	if [ ! -e "${r2_singletons}" ]; then
		touch ${r2_singletons}
	fi
	
	name=`basename ${r1} | sed 's/.atropos_final.prinseq_1.fastq//'`

# Se ainda não existe o index dentro do diretório star_index então crie o indice:
	
	if [ ! -e "${outdir}/star_index/SAindex" ]; then
		echo "Indexing genome (${refseq}) ..."
		# --genomeSAindexNbases 12 (sugestão do alinhador)
		# --sjdbOverhang 149 (sugestão do manual)	

# Criação do índice utilizando o alinhador STAR:

		STAR 	--runThreadN        ${num_threads} \  # Número de processadores a serem utilizados
			--runMode           genomeGenerate \  # Pede para gerar o índice
			--genomeFastaFiles  ${refseq} \       # Arquivo FASTA do genoma de referência
			--genomeDir         ${outdir}/star_index \   # Diretório output, que irá conter os dados gerados
			--sjdbGTFfile       ${refgtf} \       # Arquivo de anotação do genoma de referência no formato gtf
		#	--genomeSAindexNbases 12 \            # É utilizado apenas para genomas muito pequenos, escala, no caso não é necessário
			--sjdbOverhang      149 \             # Tamanho das reads -1
		 > ${outdir}/star_index/STAR.index.log.out.txt \ # Arquivo de saída com detalhes do processo
		2> ${outdir}/star_index/STAR.index.log.err.txt
	
	fi

# A função echo imprime na tela a mensagem entre "aspas"

	echo "STAR alignment PE with sample ${name}: ${r1} & ${r2} ..."
	
	# Criação do diretório para cada amostra (Conjunto de reads) dentro do diretório star_out_pe

	mkdir -p ${outdir}/star_out_pe/${name}

# Alinhamento das reads contra o genoma de referência:
	
	STAR	--runThreadN        ${num_threads} \ # Número de processadores
		--genomeDir         ${outdir}/star_index \  # Index criado na etapa anterior
		--readFilesIn       ${r1} ${r2} \  # Reads 1 e 2
		--outSAMstrandField intronMotif \  # Cria um arquivo contendo o alinhamento de splicings - É preciso para o Cufflinks
		--outFilterIntronMotifs RemoveNoncanonical \ # Exclui os alinhamento que contem junções não canonicas
		--sjdbGTFfile       ${refgtf} \  # Arquivo de anotação do genoma de referência no formato gtf
		--outFilterMultimapNmax 15 \     # Número máximo de alinhamentos multiplos permitido por read
		--outFileNamePrefix ${outdir}/star_out_pe/${name}/ \ #Prefixo dos nomes dos arquivos de saída
		--outSAMtype        BAM Unsorted \  # Tipo de arquivo de alinhamento gerado: Formato BAM, sem ordenação
		 > ${outdir}/star_out_pe/${name}/STAR.alignment_pe.log.out.txt \
		2> ${outdir}/star_out_pe/${name}/STAR.alignment_pe.log.err.txt
	
	echo "STAR alignment SE with sample ${name}: ${r1_singletons} & ${r2_singletons} ..."

# Criação dos diretórios por amostra dentro diretório star_out_se, referente aos singletons
	mkdir -p ${outdir}/star_out_se/${name}

# Mapeamento das reads singletons:

	STAR	--runThreadN        ${num_threads} \
		--genomeDir         ${outdir}/star_index \
		--readFilesIn       ${r1_singletons},${r2_singletons} \ # Reads 1 e 2 que forem singletons
		--sjdbGTFfile       ${refgtf} \
		--outSAMtype        BAM Unsorted \
		--outFilterMultimapNmax 15 \
		--outSAMstrandField intronMotif \
		--outFileNamePrefix ./$outdir/star_out_se/${name}/ \
		 > ./${outdir}/star_out_se/${name}/STAR.alignment_se.log.out.txt \
		2> ./${outdir}/star_out_se/${name}/STAR.alignment_se.log.err.txt
	
	echo "Merging STAR alignment PE & SE (${name}) ..."

# Criação do diretório star_out_final que é a junção dos resultados paired-ends e singletons:
	mkdir -p ${outdir}/star_out_final/${name}

# Combinando os resultados do alinhamento com reads paired-end e alinhamento com reads single-end (singletons)       
	 samtools merge -@ ${num_threads} -f -n  ${outdir}/star_out_final/${name}/Aligned.out.bam \
                                                ${outdir}/star_out_pe/${name}/Aligned.out.bam \
                                                ${outdir}/star_out_se/${name}/Aligned.out.bam \
	 > ${outdir}/star_out_final/${name}/samtools.merge.log.out.txt \
	2> ${outdir}/star_out_final/${name}/samtools.merge.log.err.txt

# @: Número de processadores
# -f: Permitir sobreescrever se o arquivo já existir
# -n: Ordenação por nome

	echo "Sorting STAR alignment final (${name}) ..."

# Ordenando o resultado do alinhamento por coordenadas genômicas
# - exigência para executar o cufflinks

	 samtools sort -@ ${num_threads} -o      ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \
                                                ${outdir}/star_out_final/${name}/Aligned.out.bam \
	 > ${outdir}/star_out_final/${name}/samtools.sort.log.out.txt \
	2> ${outdir}/star_out_final/${name}/samtools.sort.log.err.txt

	echo "Collecting alignment statistics (${name}) ..."

# Script para se obter as estatisticas finais:	

	SAM_nameSorted_to_uniq_count_stats.pl ${outdir}/star_out_final/${name}/Aligned.out.bam > ${outdir}/star_out_final/${name}/Aligned.stats.txt
	
	echo "Running Cufflinks (${name}) ..."

# Criação dos diretórios de cada amostra dentro do diretório cufflinks	

	mkdir -p ${outdir}/cufflinks/${name}

# Montagem dos transcritos utilizando um genoma de referência, utilizando o programa cufflinks:
	
	cufflinks --output-dir ${outdir}/cufflinks/${name} \ # Diretório output, onde ficarão os resultados
		  --num-threads ${num_threads} \    # Número de processadores
		  --GTF-guide ${refgtf} \           # Anotação de genoma de referência no formato gtf
		  --frag-bias-correct ${refseq} \   # Arquivo fasta do genoma de referência
		  --multi-read-correct \            # Para fazer uma estimação das reads que se mapeiam em multiplos locais
		  --library-type fr-unstranded \    # Tipo da biblioteca
		  --frag-len-mean 300 \             # Tamnho médio dos fragmentos
		  --frag-len-std-dev 30 \           # Desvio padrão do tamanho dos fragmentos
		  --total-hits-norm \               # Número de hits mapeados utilizando FPKM
		  --max-frag-multihits 15 \         # Número maximo de vezes que um fragmento pode se ligar a multplos locais
		  --min-isoform-fraction 0.15 \     # Porcentagem minima da abundancia da isoforma, abaixo disso são eliminadas
		  --max-intron-length 9000 \        # Tamanho máximo do intron
		  --min-intron-length 80 \          # Tamanho minimo do intron
		  --3-overhang-tolerance 300 \      # Quantidade máxima de pares de bases permitos sobrarem na extremidade 3',quando a read se alinha a referência
		  --max-bundle-frags 1000000 \      # Conjunto máximo de fragmentos em um locus
		  --max-multiread-fraction 0.55 \   # Proporção de transFrags que podem ser alinhar em multiplos locais
		  --overlap-radius 10 \             # Distancia para que os transfrags sejam unidos, se a distancia for menor que essa serão unidos
		  ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \ # Arquivo de alinhamento gerado com o STAR
		 > ${outdir}/star_out_final/${name}/cufflinks.log.out.txt \
		2> ${outdir}/star_out_final/${name}/cufflinks.log.err.txt


	echo "Running StringTie (${name}) ..."

# Criação dos diretórios para cada amostra dentro do diretório StringTie

	mkdir -p ${outdir}/stringtie/${name}

# Montagem dos transcritos utilizando genoma de referência, com o programa StringTie	

	stringtie ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \ # Alinhamento gerado pelo STAR
		-G ${refgtf} \  # Anotação do genoma de referência no formato gtf
		-o ${outdir}/stringtie/${name}/transcripts.gtf \ # Output
		-p ${num_threads} \  # Número de processadores
		-f 0.15 \ # Porcentagem minima da abundancia da isoforma, abaixo disso são eliminadas
		-a 10 \   # Junções que não apresentam reads mapedas entre ela e que não apresentem pelo menos esse valor em ambos os lados, são eliminadas
		-j 5 \    # Número minimo de reads que se alinhaem nas junções
		-c 2 \    # Mínimo de cobertura das reads
		-g 10 \   # Minimo de distancia para ser um gap em um locus.
		-M 0.45 \ # Proporção de transFrags que podem ser alinhar em multiplos locais
		-A ${outdir}/stringtie/${name}/gene_abundance.txt # output gerado: gene_abundance.txt para cada amostra


done


echo "Running cuffmerge ..."

# O comando find busca todos os transcripts.gtf gerados dentro de todos os diretórios de cada amostra.
# É gerado um arquivo txt contendo todos os arquivos encontrados

find ${outdir}/cufflinks/ -name 'transcripts.gtf' > ${outdir}/cuffmerge/assembly_GTF_list.txt

# O cuffmerge irá fundir todos os arquivos de montagem gerados para cada amostra em um único arquivo gtf.
cuffmerge -o ${outdir}/cuffmerge/ \ # Diretório de saída dos dados gerados
	--ref-gtf ${refgtf} \       # Anotação do genoma de referência em formato gtf
	--ref-sequence ${refseq} \  # Arquivo fasta do genoma de referência
	--min-isoform-fraction 0.15 \   # Porcentagem minima da abundancia da isoforma, abaixo disso são eliminadas
	--num-threads ${num_threads} \  # Número de processadores
	   ${outdir}/cuffmerge/assembly_GTF_list.txt \  # Arquivo txt gerado acima
	 > ${outdir}/cuffmerge/cuffmerge.log.out.txt \
	2> ${outdir}/cuffmerge/cuffmerge.log.err.txt

echo "Running stringtie merge ..."

# De maneira análogo o mesmo é feito para as saídas do StringTie, todos as montagens são unidas em um único arquivo:

find ${outdir}/stringtie/ -name 'transcripts.gtf' > ${outdir}/stringmerge/assembly_GTF_list.txt

stringtie --merge \
	-G ${refgtf} \  # Anotação do genoma de referência em formato gtf
	-o ${outdir}/stringmerge/merged.gtf \  # output
	-c 1 \  # Minimo de cobertura dos transcritos para serem adicionados na fusão
	-T 1 \  # Minimo de trasncritos a serem adicionados na fusão - TPM
	-f 0.15 \ # Porcentagem minima da abundancia da isoforma, abaixo disso são eliminadas
	-g 10 \   # Minimo de distancia para ser um gap em um locus.
	-i \    # Matem unidos os transcritos com retenção de introns
	${outdir}/stringmerge/assembly_GTF_list.txt

# O cuffcompare realiza uma comparação da montagem obtida com o genoma de referência, a montagem utilizada foi a saída do programa StringTie
cuffcompare	-r ${refgtf} \  # Anotação do genoma de referência em formato gtf
		-s ${refseq} \  # Arquivo fasta do genoma de referência
		-o ${outdir}/cuffcompare/stringmerge \  # Output
		${outdir}/stringmerge/merged.gtf \      # Arquivo de montagem
		 > ${outdir}/stringmerge/cuffcompare.log.out.txt \
		2> ${outdir}/stringmerge/cuffcompare.log.err.txt

######
### Using stringtie 
######

biogroup_label=()
for bamfile in `ls ${outdir}/star_out_final/*/Aligned.out.sorted.bam`; do
	name=`basename $(dirname ${bamfile})`
	echo "Running cuffquant using sample ${name} with ${outdir}/stringmerge/merged.gtf as reference ..."

# Criação dos diretórios para cada amostra, dentro do diretório cuffquant

	mkdir -p ${outdir}/cuffquant/${name}

# Quantificação dos transcritos para análise de expressão gênica diferencial

	cuffquant 	--output-dir ${outdir}/cuffquant/${name} \ # Output - Diretório de saída dos resultados
			--frag-bias-correct ${refseq} \   # Arquivo fasta do genoma de referência
			--multi-read-correct \            # Para fazer uma estimação das reads que se mapeiam em multiplos locais
			--num-threads ${num_threads} \    # Número de processadores
			--library-type fr-unstranded \    # Tipo da biblioteca
			--frag-len-mean 300 \             # Tamnho médio dos fragmentos
			--frag-len-std-dev 30 \           # Desvio padrão do tamanho dos fragmentos
			--max-bundle-frags 1000000 \      # Conjunto máximo de fragmentos em um locus
			--max-frag-multihits 15 \         # Número maximo de vezes que um fragmento pode se ligar a multiplos locais
			${outdir}/stringmerge/merged.gtf \  # Arquivo de montagem
			${bamfile} \                        # Arquivo de alinhamento no fomato BAM
		 > ${outdir}/cuffquant/${name}/cuffquant.log.out.txt \
		2> ${outdir}/cuffquant/${name}/cuffquant.log.err.txt

	groupname=`echo ${name} | sed 's/[0-9]\+$//'`
	biogroup_label=($(printf "%s\n" ${biogroup_label[@]} ${groupname} | sort -u ))

done
biogroup_files=()

echo "Running Differential Expression Analysis ..."

for label in ${biogroup_label[@]}; do
	echo -e "\tCollecting .cxb files for ${label} ..."
	group=()
	for cxbfile in `ls ${outdir}/cuffquant/${label}*/abundances.cxb`; do
		echo -e "\t\tFound ${cxbfile}"
		group=(${group[@]} "${cxbfile}")
	done
	biogroup_files=(${biogroup_files[@]} $(IFS=, ; echo "${group[*]}") )
done

echo -e "\tRunning cuffnorm & cuffdiff ..."
echo -e "\t\tLabels.: " $(IFS=, ; echo "${biogroup_label[*]}")
echo -e "\t\tFiles..: " ${biogroup_files[*]}

echo -e "\t\t\tGenerating abundance matrices (cuffnorm) ..."

# Criação do diretório cuffnorm

mkdir -p ${outdir}/cuffnorm/

# Normalização da quantificação obtida

cuffnorm 	--output-dir ${outdir}/cuffnorm \  # Diretório de saida dos dados gerados
 		--labels $(IFS=, ; echo "${biogroup_label[*]}") \  Nome de cada grupo
 		--num-threads ${num_threads} \ # Número de processadores
		--library-type fr-unstranded \  # Tipo da biblioteca
 		--library-norm-method geometric \  # Método de normalização
		--output-format simple-table \     # Formato do output
 		${outdir}/stringmerge/merged.gtf \  # Arquivo de montagem
 		${biogroup_files[*]} \              # Arquivos dos grupos separados por virgula
 	 	> ${outdir}/cuffnorm/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffnorm/cuffdiff.log.err.txt


echo -e "\t\t\tAnalysing differential expression (cuffdiff) ..."

# Criação do diretório cuffdiff

mkdir -p ${outdir}/cuffdiff/

# Análise de expressão gênica diferencial

cuffdiff 	--output-dir ${outdir}/cuffdiff \ # Diretório de saída dos dados gerados
 		--labels $(IFS=, ; echo "${biogroup_label[*]}") \  # Grupos
 		--frag-bias-correct ${refseq} \  # Arquivo fasta do genoma de referência
 		--multi-read-correct \           # Para fazer uma estimação das reads que se mapeiam em multiplos lo
 		--num-threads ${num_threads} \   # Número de processadores
 		--library-type fr-unstranded \   # tipo da biblioteca
 		--frag-len-mean 300 \            # Tamanho médio dos fragmentos
 		--frag-len-std-dev 30 \          # Desvio padrão dos fragmentos
 		--max-bundle-frags 1000000 \     # Conjunto máximo de fragmentos em um locus
                --max-frag-multihits 15 \        # Número maximo de vezes que um fragmento pode se ligar a multiplos locais
 		--total-hits-norm \              # Número de hits mapeados utilizando FPKM
 		--min-reps-for-js-test 2 \       # Número minimo de réplicas para se realizar a análise
 		--library-norm-method geometric \  # Método de normalização
 		--dispersion-method per-condition \     # Método de dispersão
 		--min-alignment-count 10 \              # Número minimo de alinhamentos em um locus para se realizar a análise
 		${outdir}/stringmerge/merged.gtf \ # Arquivo de montagem
 		${biogroup_files[*]} \
 	 	> ${outdir}/cuffdiff/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffdiff/cuffdiff.log.err.txt
 

######
### Using cufflinks/cuffmerge
######

biogroup_label=()
for bamfile in `ls ${outdir}/star_out_final/*/Aligned.out.sorted.bam`; do
	name=`basename $(dirname ${bamfile})`
	echo "Running cuffquant using sample ${name} with ${outdir}/cuffmerge/merged.gtf as reference ..."
	mkdir -p ${outdir}/cuffquant2/${name}

	cuffquant 	--output-dir ${outdir}/cuffquant2/${name} \
			--frag-bias-correct ${refseq} \
			--multi-read-correct \
			--num-threads ${num_threads} \
			--library-type fr-unstranded \
			--frag-len-mean 300 \
			--frag-len-std-dev 30 \
			--max-bundle-frags 1000000 \
			--max-frag-multihits 15 \
			${outdir}/cuffmerge/merged.gtf \
			${bamfile} \
		 > ${outdir}/cuffquant2/${name}/cuffquant.log.out.txt \
		2> ${outdir}/cuffquant2/${name}/cuffquant.log.err.txt

	groupname=`echo ${name} | sed 's/[0-9]\+$//'`
	biogroup_label=($(printf "%s\n" ${biogroup_label[@]} ${groupname} | sort -u ))

done
biogroup_files=()

echo "Running Differential Expression Analysis ..."
for label in ${biogroup_label[@]}; do
	echo -e "\tCollecting .cxb files for ${label} ..."
	group=()
	for cxbfile in `ls ${outdir}/cuffquant2/${label}*/abundances.cxb`; do
		echo -e "\t\tFound ${cxbfile}"
		group=(${group[@]} "${cxbfile}")
	done
	biogroup_files=(${biogroup_files[@]} $(IFS=, ; echo "${group[*]}") )
done

echo -e "\tRunning cuffnorm & cuffdiff ..."
echo -e "\t\tLabels.: " $(IFS=, ; echo "${biogroup_label[*]}")
echo -e "\t\tFiles..: " ${biogroup_files[*]}

echo -e "\t\t\tGenerating abundance matrices (cuffnorm) ..."

mkdir -p ${outdir}/cuffnorm2/

cuffnorm 	--output-dir ${outdir}/cuffnorm2 \
 		--labels $(IFS=, ; echo "${biogroup_label[*]}") \
 		--num-threads ${num_threads} \
		--library-type fr-unstranded \
 		--library-norm-method geometric \
		--output-format simple-table \
 		${outdir}/cuffmerge/merged.gtf \
 		${biogroup_files[*]} \
 	 	> ${outdir}/cuffnorm2/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffnorm2/cuffdiff.log.err.txt


echo -e "\t\t\tAnalysing differential expression (cuffdiff) ..."

mkdir -p ${outdir}/cuffdiff2/

cuffdiff 	--output-dir ${outdir}/cuffdiff2 \
 		--labels $(IFS=, ; echo "${biogroup_label[*]}") \
 		--frag-bias-correct ${refseq} \
 		--multi-read-correct \
 		--num-threads ${num_threads} \
 		--library-type fr-unstranded \
 		--frag-len-mean 300 \
 		--frag-len-std-dev 50 \
 		--max-bundle-frags 9999999 \
 		--max-frag-multihits 20 \
 		--total-hits-norm \
 		--min-reps-for-js-test 2 \
 		--library-norm-method geometric \
 		--dispersion-method per-condition \
 		--min-alignment-count 10 \
 		${outdir}/cuffmerge/merged.gtf \
 		${biogroup_files[*]} \
 	 	> ${outdir}/cuffdiff2/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffdiff2/cuffdiff.log.err.txt
 


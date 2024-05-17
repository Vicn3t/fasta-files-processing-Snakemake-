# Snakefile assignment Vincent Spinelli

# command to execute : sbatch mysbatch.slurm

# Configuring the constants
refGenome = "/lustre1/project/stg_00079/teaching/hg38_21/chr21.fa"
raw_sample_names, = glob_wildcards("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz")

# filtering the list of sample names by keeping names with a 7 on its 6th character:
sample_names = [name for name in raw_sample_names if name[6] == '7']

#Defining all the files we want to produce
rule all:
    input:
        fastq_unzipped=expand("00.unzipped_fastq/{sample}.fastq", sample=sample_names),

        fastqc_zip=expand("01.fastqc/{sample}_fastqc.zip", sample=sample_names),
        fastqc_report = expand("01.fastqc/{sample}_fastqc.html", sample=sample_names),
        summary_data=expand("01.fastqc/{sample}_fastqc/fastqc_data.txt", sample=sample_names),

        bam_file = expand("02.bam/{sample}.bam", sample=sample_names),
        bai_file = expand("02.bam/{sample}.bam.bai", sample=sample_names),

        raw_vcf = "03.raw_snps/all_samples_raw.vcf",

        cleaned_vcf = "04.cleaned_snps/all_samples_cleaned.vcf",

        annotated_vcf = "05.annotated_snps/all_samples_annotated.vcf",

        variant_quality_histogram_raw = "03.raw_snps/variant_quality_histogram_raw.png",
        variant_quality_histogram_cleaned = "04.cleaned_snps/variant_quality_histogram_cleaned.png",

        genes_vcf = "genes.vcf",

        heatmap = "05.annotated_snps/snp_gene_heatmap.png",


# Retrieving and unzipping the fastq files from the 1000genomes folder
rule unzipped_fastq:
    input:
        "/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz"
    output:
        unzipped_fastq = "00.unzipped_fastq/{sample}.fastq"
    shell:
        """
        mkdir -p 00.unzipped_fastq

        # Checking the file being unzipped
        echo "Unzipping of {output}"

        gunzip -c {input} > {output}

        # Checking if the file has been generated
        if [ ! -f {output.unzipped_fastq} ]; then
        echo "Error: Unzipping for {wildcards.sample} failed."
        exit 1
        fi

        # Keeping track on the file that has finished unzipping
        echo "Unzipping of {output} completed"
        """


# Fastqc analysis of the unzipped fastq files
rule fastqc:
    input:
        fastq="00.unzipped_fastq/{sample}.fastq"
    output:
        fastqc_zip="01.fastqc/{sample}_fastqc.zip",
        fastqc_report="01.fastqc/{sample}_fastqc.html",
        summary_data="01.fastqc/{sample}_fastqc/fastqc_data.txt",
    shell:
        """
        mkdir -p 01.fastqc

        echo "Fastqc analysis of {input.fastq} "

        fastqc {input.fastq} -o 01.fastqc --extract
        
        if [ ! -f {output.fastqc_zip} ]; then
        echo "Error: FastQC analysis for {wildcards.sample} failed."
        exit 1
        fi

        echo "FastQC analysis of {input.fastq} completed"
        """


# Mapping the fastq files to the reference genome with their bwa index
rule bam_bai:
    input:
        fastq="00.unzipped_fastq/{sample}.fastq",
    output:
        bam = "02.bam/{sample}.bam",
        bai = "02.bam/{sample}.bam.bai",
    params:
        DB = refGenome,
    shell:
        """
        mkdir -p 02.bam

        echo "Mapping of {input.fastq} to the reference genome"

        bwa mem {params.DB} {input.fastq} | samtools sort -o {output.bam} && samtools index {output.bam}

        echo "Mapping of {input.fastq} to the reference genome completed"
        """


# Variant calling on the bam files 
rule variant_calling:
    input:
        bams=expand("02.bam/{sample}.bam", sample = sample_names),
    output:
        vcf="03.raw_snps/all_samples_raw.vcf",
    params:
        DB=refGenome,
    shell:
        """
        mkdir -p 03.raw_snps

        echo "variants calling on the bam files"

        bcftools mpileup -Ou -f {params.DB} {input.bams} | \
        bcftools call -mv -Ov -o {output.vcf}

        if [ ! -f 03.raw_snps/all_samples_raw.vcf ]; then
        echo "Error: variant calling step failed."
        exit 1
        fi

        echo "variants calling on the bam files completed"
      """
    

# Variant cleaning the vcf file with filtering and normalization 
rule variant_cleanup:
    input:
        vcf="03.raw_snps/all_samples_raw.vcf"
    output:
        vcf="04.cleaned_snps/all_samples_cleaned.vcf"
    params:
        db=refGenome,
    shell:
        """
        mkdir -p 04.cleaned_snps

        echo "cleaning {input.vcf}"

        ( cat {input.vcf} \
           | vt decompose - \
           | vt normalize -n -r {params.db} - \
           | vt uniq - \
           | vt view -f "QUAL>20" -h - \
           > {output.vcf} )

        if [ ! -f 04.cleaned_snps/all_samples_cleaned.vcf ]; then
        echo "Error: Cleaning step failed."
        exit 1
        fi

        echo "cleaning {input.vcf} completed"
        """


# Annotating the vcf file using snpeff
rule snpeff_annotation:
    input:
        vcf="04.cleaned_snps/all_samples_cleaned.vcf",
    params:
        # importing the snpeff components from the config file
        snpeff_db_folder = "/staging/leuven/stg_00079/teaching/snpeff_db",
        snpeff_jar = "/lustre1/project/stg_00079/teaching/I0U19a_conda_2024/share/snpeff-5.2-0/snpEff.jar",
        snpeff_genome = "hg38",
    output:
        annotated_vcf = "05.annotated_snps/all_samples_annotated.vcf",
    shell:
        """
        mkdir -p 05.annotated_snps

        echo "Annotating the vcf file"

        java -Xmx4096m -jar \
            {params.snpeff_jar} eff {params.snpeff_genome} \
            -dataDir {params.snpeff_db_folder} \
            {input.vcf} > {output.annotated_vcf}
            mv snpEff_genes.txt snpEff_summary.html 05.annotated_snps

        if [ ! -f 05.annotated_snps/all_samples_annotated.vcf ]; then
        echo "Error: Annotation step failed."
        exit 1
        fi

        echo "Annotating the vcf completed"

        """


# Extracting SNPs associated with APP/SOD1/DYRK1A from the VCF file and save it to the root folder
rule extract_gene_specific_snps:
    input:
        vcf = "05.annotated_snps/all_samples_annotated.vcf"
    output:
        genes_vcf = "genes.vcf"
    shell:
        """

        mkdir -p .

        echo "Extracting SNPs associated with specific genes from {input.vcf}"

        # Extracting the header
        grep '^#' {input.vcf} > {output.genes_vcf}

        # Extracting the SNPs associated with APP, SOD1, and DYRK1A
        bcftools view {input.vcf} | \
        awk '/ANN=.*(APP|SOD1|DYRK1A)/' >> {output.genes_vcf}

        if [ ! -f genes.vcf ]; then
        echo "Error: genes specific selection step failed."
        exit 1
        fi

        echo "Extraction of SNPs associated with specific genes completed"
        """

# Putting into the report a histogram of the variants quality scores before cleaning
rule variant_quality_histogram_raw:
    input:
        vcf = "03.raw_snps/all_samples_raw.vcf"
    output:
        variant_quality_histogram_raw = report("03.raw_snps/variant_quality_histogram_raw.png",
                           category='Variant Quality',
                           subcategory='Quality Histogram',
                           labels={"Description": "Raw Data Quality Score Distribution"})
    run:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        import os

        # Create the directory if it does not exist
        os.makedirs("03.raw_snps", exist_ok=True)

        # Extract the quality scores from the VCF file
        quality_scores = []
        with open(input.vcf, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    qual = float(fields[5])
                    quality_scores.append(qual)

        data = pd.DataFrame(quality_scores, columns=['Quality'])

        # Plot the histogram
        plt.figure(figsize=(6, 5))
        sns.histplot(data['Quality'], kde=True)
        plt.title("Raw Data Distribution of Quality Scores")
        plt.xlabel("Quality")
        plt.ylabel("Frequency")
        plt.xlim(0, 500)  
        plt.tight_layout()  
        plt.savefig(output.variant_quality_histogram_raw)  

# Putting into the report a histogram of the variants quality scores after cleaning
rule variant_quality_histogram_cleaned:
    input:
        vcf = "04.cleaned_snps/all_samples_cleaned.vcf"
    output:
        variant_quality_histogram_cleaned = report("04.cleaned_snps/variant_quality_histogram_cleaned.png",
                           category='Variant Quality',
                           subcategory='Quality Histogram',
                           labels={"Description": "Cleaned Data Quality Score Distribution"})
    run:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        import os

        # Create the directory if it does not exist
        os.makedirs("04.cleaned_snps", exist_ok=True)

        # Extract the quality scores from the VCF file
        quality_scores = []
        with open(input.vcf, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    qual = float(fields[5])
                    quality_scores.append(qual)


        # Create a DataFrame with the quality scores
        data = pd.DataFrame(quality_scores, columns=['Quality'])

        # Plot the histogram
        plt.figure(figsize=(6, 5))
        sns.histplot(data['Quality'], kde=True)
        plt.title("Cleaned Data Distribution of Quality Scores")
        plt.xlabel("Quality")
        plt.ylabel("Frequency")
        plt.xlim(0, 500)  
        plt.tight_layout()  
        plt.savefig(output.variant_quality_histogram_cleaned)  

# Generating a heatmap of SNPs associated with APP, SOD1, and DYRK1A per individual
rule gene_snp_heatmap:
    input:
        genes_vcf = "genes.vcf",
    output:
        heatmap = report("05.annotated_snps/snp_gene_heatmap.png",
                         category='Gene Heatmap',
                         subcategory='SNP Associations',
                         labels={"Description": "SNPs associated with specific genes per individual."}),
    run:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        import os
        import re

        os.makedirs("05.annotated_snps", exist_ok=True)

        gene_snp_counts = {"APP": {}, "SOD1": {}, "DYRK1A": {}}

        # Extract the gene information and count SNPs for each individual
        with open(input.genes_vcf, 'r') as f:
            header = None
            for line in f:
                if line.startswith('#CHROM'):
                    # Extract individual identifiers from the header
                    header = [re.sub(r'^02\.bam/|\.GRCh38DH\.exome\.chr21\.bam$', '', sample) for sample in line.strip().split('\t')[9:]]
                elif not line.startswith('#'):
                    fields = line.strip().split('\t')
                    gene_annotation = re.search(r'ANN=.*(APP|SOD1|DYRK1A)', fields[7])
                    if gene_annotation:
                        gene = gene_annotation.group(1)
                        # Count SNPs for each individual associated with this gene
                        for i, sample in enumerate(header):
                            genotype_info = fields[9 + i].split(':')[0] 
                            if genotype_info not in ('./.', '0/0'):
                                if sample not in gene_snp_counts[gene]:
                                    gene_snp_counts[gene][sample] = 0
                                gene_snp_counts[gene][sample] += 1

        df = pd.DataFrame(gene_snp_counts).fillna(0)

        plt.figure(figsize=(12, 8))
        sns.heatmap(df, annot=True, fmt=".0f", cmap='coolwarm', cbar=True)
        plt.title("Heatmap of SNPs Associated with Genes per Individual")
        plt.xlabel("Genes")
        plt.ylabel("Individuals")


        plt.yticks(rotation=0)

        plt.tight_layout()

        plt.savefig(output.heatmap)

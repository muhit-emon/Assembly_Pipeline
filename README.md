# Assembly_Pipeline
Assembly pipeline of CIWARS for computing resistome risk and identifying ARG like sequencs and their mobility
# Requirements
<ol>
  <li>Linux operating system</li>
  <li>conda</li>
</ol>

# Installation
<pre>
git clone https://github.com/muhit-emon/Assembly_Pipeline.git
cd Assembly_Pipeline
bash install.sh
conda env create -f environment.yml
</pre>
# conda environment activation
After installation, a conda environment named <b>assembly_pipeline</b> will be created.<br>
To activate the environment, run the following command <br>
<pre>
conda activate assembly_pipeline
</pre>
# Download the compressed Blast Database file from Zenodo (25 GB) and uncompress it
Go inside <b>Assembly_Pipeline</b> directory
<pre>
wget https://zenodo.org/records/10471551/files/BlastDB.tar.gz
tar -zxvf BlastDB.tar.gz
</pre>

# Download the compressed DeepARG-DB and mobileOG databases (DB.tar.gz) file from one drive, put it inside the "Assembly_Pipeline" directory and uncompress it
Go to <a href="https://virginiatech-my.sharepoint.com/:u:/g/personal/muhitemon_vt_edu/EQjIpLhOmMVFotPRDiK5Id0BGQFVSVVbyEWOmRWiz-rYUA">One Drive</a> and download DB.tar.gz. Put it inside the <b>Assembly_Pipeline</b> directory
<pre>
tar -zxvf DB.tar.gz
</pre>

# (1) Usage on metagenomic paired-end short read data
Go inside <b>meta-VF-AMR/scripts</b> directory. <br> <br>
<b>To run metaVF-AMR on metagenomic paired-end short read data (<span> &#42; </span>.fastq/<span> &#42; </span>.fq/<span> &#42; </span>.fastq.gz/<span> &#42; </span>.fq.gz), use the following command</b> <br>
<pre>
nextflow run pe_pipeline.nf --R1 &ltabsolute/path/to/forward/read/file&gt --R2 &ltabsolute/path/to/reverse/read/file&gt --out_fname &ltprefix of output file name&gt
rm -r work
</pre>
The command line options for this script (<b>pe_pipeline.nf</b>) are: <br><br>
<b>--R1</b>: The absolute path of the fastq file containing forward read sequences <br>
<b>--R2</b>: The absolute path of the fastq file containing reverse read sequences <br>
<b>--out_fname</b>: The prefix of the output file name <br><br>

With <b>--out_fname demo</b>, Four fasta files named <b>demo_Reconstructed_ARG.fasta</b>, <b>demo_Reconstructed_MRG.fasta</b>, <b>demo_Reconstructed_BRG.fasta</b>, and <b>demo_Reconstructed_VF.fasta</b> will be generated inside <b>meta-VF-AMR/TMP</b> directory. <br><br>

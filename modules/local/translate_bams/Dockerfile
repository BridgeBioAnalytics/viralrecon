
FROM nfcore/base:2.1
LABEL \
    author="Olga Botvinnik" \
    description="Image for capsid evolution translate bam" \
    maintainer="olga.botvinnik@bridgebio.com"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-translate-bam/bin:$PATH


# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-translate-bam > nf-core-translate-bam.yml

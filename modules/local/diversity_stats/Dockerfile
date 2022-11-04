
FROM nfcore/base:2.1
LABEL \
    author="Olga Botvinnik" \
    description="Image for capsid evolution compute diversity" \
    maintainer="olga.botvinnik@bridgebio.com"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-compute-diversity/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-compute-diversity > nf-core-compute-diversity.yml

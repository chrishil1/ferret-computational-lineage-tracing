# Get base image
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# environment
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH $PATH:/opt/conda/bin
# ENV RENV_VERSION 0.16.0

# get conda installed
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean --all -f -y && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    apt-get update

RUN /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

COPY environment.yaml .
RUN /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    conda env update -n base -f environment.yaml && \
    conda clean -afy

# Install R packages with renv
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"

WORKDIR /home/rstudio
COPY R_dependencies.R R_dependencies.R
COPY renv.lock renv.lock

# Install Fix Bioconductor + spatial ecosystem
RUN R -e "\
install.packages(c('spatstat.utils','spatstat.data','spatstat.geom','ggplot2'), \
repos='https://cloud.r-project.org'); \
BiocManager::install(c('ggtree','treeio','clusterProfiler','enrichplot'), \
ask=FALSE, update=TRUE)"

# Install Seurat core
RUN R -e "install.packages('Seurat', repos='https://cloud.r-project.org')"
RUN R -e "renv::restore(prompt=FALSE, exclude=c('Seurat','SeuratData','SeuratDisk','ggtree','clusterProfiler','enrichplot','scAnalysisR'))"

# Install Seurat ecosystem packages last
RUN R -e "\
install.packages('remotes', repos='https://cloud.r-project.org'); \
remotes::install_github('satijalab/seurat-data'); \
remotes::install_github('mojaveazure/seurat-disk'); \
remotes::install_github('CUAnschutzBDC/scAnalysisR'); \
stopifnot(all(c('SeuratData','SeuratDisk','scAnalysisR') %in% installed.packages()[,1]))"

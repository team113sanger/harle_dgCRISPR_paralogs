FROM rocker/r-base:4.2.3

LABEL maintainer="vo1@sanger.ac.uk" \
      version="1.0.0" \
      description="Paralog dgCRISPR1"

MAINTAINER  Victoria Offord <vo1@sanger.ac.uk>

USER root

# Locale
ENV LC_ALL C
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

# Environment variables
ENV VIRTUAL_ENV='/opt/venv'
ENV VER_MAGECK_MAJOR="0.5"
ENV VER_MAGECK_MINOR="9.3"
ENV BAGEL_INST_DIR="/opt/bagel"
ENV VER_BAGEL="f9eedca"

# Prevent interactive options
ENV DEBIAN_FRONTEND=noninteractive

# Update packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \ 
    libv8-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    curl \
    software-properties-common \
    gcc \
    git \
    cmake \
    && rm -rf /var/lib/apt/lists/* 

# Install R library dependencies
RUN R -q -e "install.packages('curl')" \
    && strip /usr/local/lib/R/site-library/*/libs/*.so
RUN R -q -e "install.packages('optparse')"
RUN R -q -e "install.packages('stringi')"
RUN R -q -e "install.packages('tidyverse')"
RUN R -q -e "install.packages('reshape2')"
RUN R -q -e "install.packages('vroom')"
RUN R -q -e "install.packages('data.table')"
RUN R -q -e "install.packages('ggpubr')"
RUN R -q -e "install.packages('ggsci')"
RUN R -q -e "install.packages('ggrepel')"
RUN R -q -e "install.packages('patchwork')"
RUN R -q -e "install.packages('GGally')"
RUN R -q -e "install.packages('ggridges')"
RUN R -q -e "install.packages('pROC')"
RUN R -q -e "install.packages('gridExtra')"
RUN R -q -e "install.packages('VennDiagram')"

# Install Python 3.8 as 3.10 has issues with dependencies (numpy)
RUN wget https://www.python.org/ftp/python/3.8.0/Python-3.8.0.tgz \
    && tar -xf Python-3.8.0.tgz \
    && cd Python-3.8.0 \
    && ./configure \
    && make -j 12 \
    && make altinstall \
    && cd .. \
    && rm -rf Python-3.8.0*

# Copy Python requirements file
COPY requirements.txt requirements.txt

# Create Python 3.8 virtual environment
RUN python3.8 -m venv "${VIRTUAL_ENV}" && chmod +x "${VIRTUAL_ENV}"/bin/activate && "${VIRTUAL_ENV}"/bin/activate && "${VIRTUAL_ENV}"/bin/pip3 install --no-cache-dir -r requirements.txt 

# Install MAGeCK
RUN curl -sSL --retry 10 -o mageck.tar.gz https://downloads.sourceforge.net/project/mageck/${VER_MAGECK_MAJOR}/mageck-${VER_MAGECK_MAJOR}.${VER_MAGECK_MINOR}.tar.gz \
    && mkdir mageck \
    && tar --strip-components 1 -C mageck -xzf mageck.tar.gz \
    && cd mageck \
    && "${VIRTUAL_ENV}"/bin/activate \
    && "${VIRTUAL_ENV}"/bin/python3 setup.py install \
    && rm -rf mageck.* mageck/*

# Install BAGEL2
RUN mkdir $BAGEL_INST_DIR \
    && git clone https://github.com/hart-lab/bagel.git ${BAGEL_INST_DIR} \
    && cd ${BAGEL_INST_DIR} \
    && git reset --hard ${VER_BAGEL}

# Add Python virtual environment and BAGEL install to PATH
ENV PATH="${VIRTUAL_ENV}/bin:/opt/bagel:$PATH"

# Set display
ENV DISPLAY=:0

# USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]

FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update
RUN apt-get update && apt-get install -y automake build-essential bzip2 wget git default-jre unzip

ENV PATH="/usr/local/bin:${PATH}"

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
 bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /usr/local/ && \
 rm Miniconda3-latest-Linux-x86_64.sh

RUN conda install -y -c conda-forge hdf5 pytables pypandoc biopython networkx numpy pandas scipy scikit-learn psutil

RUN pip install setuptools-markdown

RUN git clone https://bolduc@bitbucket.org/MAVERICLab/vcontact2.git && cd vcontact2 && \
pip install --no-dependencies .

RUN wget http://micans.org/mcl/src/mcl-latest.tar.gz && \
 tar xf mcl-latest.tar.gz && cd mcl-14-137 && ./configure --prefix /usr/local/ && make install && \
 rm -rf /mcl-latest.tar.gz

RUN wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar -P /usr/local/bin

RUN echo -e '#!/bin/bash\njava -jar /usr/local/bin/cluster_one-1.0.jar "$@"\n' > /usr/local/bin/cluster_one-1.0.sh && \
chmod +x $PATH/cluster_one-1.0.sh

RUN wget --no-verbose ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz && \
 tar xf ncbi-blast-2.6.0+-x64-linux.tar.gz && cd ncbi-blast-2.6.0+ && \
 cp bin/* $PATH && cd / && rm -rf ncbi-blast-2.6.0+*

RUN wget --no-verbose http://github.com/bbuchfink/diamond/releases/download/v0.9.10/diamond-linux64.tar.gz && \
 tar xf diamond-linux64.tar.gz && cp diamond $PATH && rm -rf diamond-linux64.tar.gz

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]

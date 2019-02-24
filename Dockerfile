FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update && apt-get install -y automake build-essential bzip2 wget git default-jre unzip

RUN add-apt-repository ppa:george-edison55/cmake-3.x && apt-get update && apt-get install -y cmake

ENV PATH="/usr/local/bin:${PATH}"

RUN pip install pandas

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
 bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /usr/local/ && \
 rm Miniconda3-latest-Linux-x86_64.sh

RUN conda install -y -c conda-forge hdf5 pytables pypandoc biopython networkx numpy pandas scipy scikit-learn psutil
RUN conda install -y -c bioconda mcl blast diamond

RUN conda clean --yes --tarballs --packages --source-cache

RUN pip install setuptools-markdown

RUN git clone https://bolduc@bitbucket.org/MAVERICLab/vcontact2.git && cd vcontact2 && \
 pip install .

RUN cp vcontact2/vcontact/data/ViralRefSeq-prokaryotes-v8?.* /usr/local/lib/python3.7/site-packages/vcontact/data/

RUN wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar -P /usr/local/bin

RUN echo -e '#!/bin/bash\njava -jar /usr/local/bin/cluster_one-1.0.jar "$@"\n' > /usr/local/bin/cluster_one-1.0.sh && \
chmod +x /usr/local/bin/cluster_one-1.0.sh

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]

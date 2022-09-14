FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

ENV PATH=/miniconda/bin:/usr/local/bin:${PATH}

# Install system dependencies
RUN apt-get update && apt-get install -y automake build-essential bzip2 wget git default-jre unzip ca-certificates

RUN conda update -y -n base conda

RUN conda install -n base -c conda-forge 'python=3.7.*' conda-build

RUN conda install -y -c bioconda -c conda-forge "vcontact2==0.9.19" hdf5 pytables pypandoc biopython networkx "numpy<1.19.5" \
    "pandas==0.25.3" scipy scikit-learn psutil pyparsing mcl blast prodigal

RUN conda clean -y --all

RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.36/diamond-linux64.tar.gz && tar xzf diamond-linux64.tar.gz
RUN chmod +x diamond && cp diamond /usr/local/bin/

RUN pip install jsonrpcbase nose jinja2 setuptools-markdown configparser #pyparsing pandas

RUN wget http://paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar && mv cluster_one-1.0.jar /usr/local/bin/
RUN chmod +x /usr/local/bin/cluster_one-1.0.jar

# Clean stuff up to reduce disk space
RUN apt-get clean && conda build purge-all

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]

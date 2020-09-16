FROM continuumio/miniconda3:4.8.2

RUN conda update -n base -c defaults conda
COPY environment.yml /tmp/environment.yml
RUN conda env create -n sarscov2phylo -f /tmp/environment.yml

RUN mkdir -p /work/sarscov2phylo
ADD . /work/sarscov2phylo
WORKDIR /work/sarscov2phylo

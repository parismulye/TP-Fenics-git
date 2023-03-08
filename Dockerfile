FROM jupyter/scipy-notebook:cf6258237ff9

USER root
RUN apt-get update -y 
RUN apt-get install -y git cmake libblkid-dev e2fslibs-dev libboost-all-dev libaudit-dev libeigen3-dev libmpc-dev

RUN python3 -m pip install --no-cache-dir notebook jupyterlab
RUN pip install --no-cache-dir jupyterhub
RUN pip install fenics-ffc==2019.1.0 --upgrade

RUN git clone --branch=2019.1.0 https://bitbucket.org/fenics-project/dolfin
RUN git clone --branch=2019.1.0 https://bitbucket.org/fenics-project/mshr
RUN mkdir dolfin/build && cd dolfin/build && cmake .. && make install && cd ../..
RUN mkdir mshr/build   && cd mshr/build   && cmake .. && make install && cd ../..


RUN pip install "pybind11[global]"

RUN cd dolfin/python && pip install . && cd ../..
RUN cd mshr/python   && pip install . && cd ../..

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}


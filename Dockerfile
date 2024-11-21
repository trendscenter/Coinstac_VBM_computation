FROM python:3.8

# Initiation of system
RUN export MCR_CACHE_VERBOSE=true
RUN apt-get update -y \
 && apt-get install -y wget unzip libxext-dev libxt-dev libxmu-dev libglu1-mesa-dev libxrandr-dev build-essential gcc

WORKDIR /computation
#-------------------------------------------------
# Install Miniconda, and set up Python environment
#-------------------------------------------------
# Install required libraries
#RUN sed -i 's|^mirrorlist=|#mirrorlist=|g' /etc/yum.repos.d/CentOS-Base.repo && \
#    sed -i 's|^#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-Base.repo
#RUN yum install -y -q gcc libXext.x86_64 libXt.x86_64 unzip zip wget \
#    && yum clean packages \
#    && rm -rf /var/cache/yum/* /tmp/* /var/tmp/*
#ENV PATH=/opt/miniconda/envs/default/bin:$PATH
#RUN curl -ssL -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
#    && bash miniconda.sh -b -p /opt/miniconda \
#    && rm -f miniconda.sh \
#    && /opt/miniconda/bin/conda update -n base conda \
#    && /opt/miniconda/bin/conda config --add channels conda-forge \
#    && /opt/miniconda/bin/conda create -y -q -n default python=3.11 \
#    && rm -rf /opt/miniconda/[!envs]*


# Install MATLAB MCR in /opt/mcr/
ENV MATLAB_VERSION R2022b
ENV MCR_VERSION R2022b
#ENV MCR_VERSION v97
ENV MCR_UPDATE 9
RUN mkdir /opt/mcr_install \
 && mkdir /opt/mcr \
 && wget --progress=bar:force -P /opt/mcr_install https://ssd.mathworks.com/supportfiles/downloads/${MATLAB_VERSION}/Release/${MCR_UPDATE}/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_${MATLAB_VERSION}_Update_${MCR_UPDATE}_glnxa64.zip \
 && unzip -q /opt/mcr_install/MATLAB_Runtime_${MATLAB_VERSION}_Update_${MCR_UPDATE}_glnxa64.zip -d /opt/mcr_install \
 && /opt/mcr_install/install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent \
 && rm -rf /opt/mcr_install /tmp/*

# Install SPM Standalone in /opt/spm12/
ENV SPM_VERSION 12
ENV SPM_REVISION r7771
ENV LD_LIBRARY_PATH /opt/mcr/${MCR_VERSION}/runtime/glnxa64:/opt/mcr/${MCR_VERSION}/bin/glnxa64:/opt/mcr/${MCR_VERSION}/sys/os/glnxa64:/opt/mcr/${MCR_VERSION}/sys/opengl/lib/glnxa64:/opt/mcr/${MCR_VERSION}/extern/bin/glnxa64
ENV MCR_INHIBIT_CTF_LOCK 1
ENV SPM_HTML_BROWSER 0
# Running SPM once with "function exit" tests the succesfull installation *and*
# extracts the ctf archive which is necessary if singularity is going to be
# used later on, because singularity containers are read-only.
# Also, set +x on the entrypoint for non-root container invocations
RUN wget --no-check-certificate --progress=bar:force -P /opt https://www.fil.ion.ucl.ac.uk/spm/download/restricted/utopia/spm12/spm${SPM_VERSION}_${SPM_REVISION}_Linux_${MATLAB_VERSION}.zip \
 && unzip -q /opt/spm${SPM_VERSION}_${SPM_REVISION}_Linux_${MATLAB_VERSION}.zip -d /opt \
 && rm -f /opt/spm${SPM_VERSION}_${SPM_REVISION}_Linux_${MATLAB_VERSION}.zip \
 && sed -i '33,35d' /opt/spm12/run_spm12.sh \
 && /opt/spm${SPM_VERSION}/spm${SPM_VERSION} function exit \
 && chmod +x /opt/spm${SPM_VERSION}/spm${SPM_VERSION}

#Unzip /opt/spm12/spm12.ctf file 
WORKDIR /opt/spm12
RUN unzip spm12.ctf \
 && chmod -R 777 /opt/spm12/fsroot 

# Set the working directory
WORKDIR /computation

# Copy the current directory contents into the container
COPY requirements.txt /computation

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy the current directory contents into the container
COPY . /computation

# allow all inputs on vbm interface
COPY preprocess.py /usr/local/lib/python3.8/site-packages/nipype/interfaces/spm/
#COPY preprocess.py /opt/miniconda/envs/default/lib/python3.11/site-packages/nipype/interfaces/spm/

CMD ["python", "/computation/entry.py"]

FROM centos:7

WORKDIR /computation
ADD . /computation

#----------------------------
# Install common dependencies
#----------------------------
RUN yum install -y -q bzip2 ca-certificates curl unzip \
    && yum clean packages \
    && rm -rf /var/cache/yum/* /tmp/* /var/tmp/*

#-------------------------------------------------
# Install Miniconda, and set up Python environment
#-------------------------------------------------
ENV PATH=/opt/miniconda/envs/default/bin:$PATH
RUN curl -ssL -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash miniconda.sh -b -p /opt/miniconda \
    && rm -f miniconda.sh \
    && /opt/miniconda/bin/conda update -n base conda \
    && /opt/miniconda/bin/conda config --add channels conda-forge \
    && /opt/miniconda/bin/conda create -y -q -n default python=3.5.1 \
    	traits pandas \
#    && conda clean -y --all \
    && pip install -U -q --no-cache-dir pip \
    && pip install -q --no-cache-dir \
    	nipype \
    && rm -rf /opt/miniconda/[!envs]*

#----------------------
# Install MCR and SPM12
#----------------------
# Install required libraries
RUN yum install -y -q libXext.x86_64 libXt.x86_64 \
    && yum clean packages \
    && rm -rf /var/cache/yum/* /tmp/* /var/tmp/*

# Install MATLAB Compiler Runtime
WORKDIR /opt
#RUN curl -sSL -o mcr.zip https://www.mathworks.com/supportfiles/downloads/R2017a/deployment_files/R2017a/installers/glnxa64/MCR_R2017a_glnxa64_installer.zip \
RUN unzip -q /computation/softwares/mcr_install/MCR_R2017a_glnxa64_installer.zip -d mcrtmp \
    && mcrtmp/install -destinationFolder /opt/mcr -mode silent -agreeToLicense yes \
    && rm -rf mcrtmp /computation/softwares/mcr_install/MCR_R2017a_glnxa64_installer.zip /tmp/*

# Install standalone SPM
WORKDIR /opt
#RUN curl -sSL -o spm.zip http://www.fil.ion.ucl.ac.uk/spm/download/restricted/utopia/dev/spm12_latest_Linux_R2017a.zip \
RUN unzip -q /computation/softwares/spm_install/spm12_latest_Linux_R2017a.zip \
    && rm -rf /computation/softwares/spm_install/spm12_latest_Linux_R2017a.zip \
    && unzip /opt/spm12/spm12.ctf -d /opt/spm12/
ENV MATLABCMD=/opt/mcr/v*/toolbox/matlab \
    SPMMCRCMD="/opt/spm*/run_spm*.sh /opt/mcr/v*/ script" \
    FORCE_SPMMCR=1 \
    LD_LIBRARY_PATH=/opt/mcr/v*/runtime/glnxa64:/opt/mcr/v*/bin/glnxa64:/opt/mcr/v*/sys/os/glnxa64:$LD_LIBRARY_PATH
    
# Install Bids Validator using npm
RUN curl -sL https://rpm.nodesource.com/setup_6.x | bash - \
    && yum install -y nodejs gcc \
    && npm install -g bids-validator 

#Install other python packages
RUN pip install med2image pybids coverage
RUN pip install --no-cache-dir -r /computation/requirements.txt

#Install bottle
ADD server/requirements.txt /server
ADD server/server.py /server
RUN pip install cython future pillow
RUN pip install --no-cache-dir -r /computation/server/requirements.txt

#Remove user warning from dicom init file
RUN sed -i '53d' /opt/miniconda/envs/default/lib/python3.5/site-packages/dicom/__init__.py

CMD ["python", "-u", "/computation/server/server.py"]

# Window in Wall



### 1. Setting up the Anaconda environment with COMPAS

Execute the commands below in Anaconda Prompt:
	
    (base) conda config --add channels conda-forge
    (base) conda create -n wiw compas_fab --yes
    (base) conda activate wiw
    
#### Verify Installation
    (wiw) pip show compas_fab

####
    Name: compas-fab
    Version: 0.xx
    Summary: Robotic fabrication package for the COMPAS Framework
    ...

#### Install on Rhino

    (wiw) python -m compas_rhino.install -v 7.0


### 2. Installation of Dependencies

    (wiw) conda install git

#### AM Information Model
    
    (wiw) python -m pip install git+https://github.com/augmentedfabricationlab/am_information_model@master#egg=am_information_model
    (wiw) python -m compas_rhino.install -p am_information_model -v 7.0

#### UR Fabrication Control
    
    (wiw) python -m pip install git+https://github.com/augmentedfabricationlab/ur_fabrication_control@master#egg=ur_fabrication_control
    (wiw) python -m compas_rhino.install -p ur_fabrication_control -v 7.0


### 3. Cloning the Course Repository

Create a workspace directory:

C:\Users\YOUR_USERNAME\workspace

Then open Github Desktop and clone the [window_in_wall](https://github.com/augmentedfabricationlab/window_in_wall) repository into you workspace folder:


### Credits

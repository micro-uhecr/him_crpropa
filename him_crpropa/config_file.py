"""A configuration file for the hadronic module

02.12.2021 - Leonel Morejon
"""

mtag = 'SIBYLL23D'   
# mtag = 'EPOSLHC'     # 1000 times slower
# mtag = 'PHOJET112'   
# mtag = 'URQMD34'     # 100 times slower
# mtag = 'PYTHIA8'     # weird behavior
# mtag = 'QGSJET01C'    

impy_directory = '/home/leonel/GitProjects/impy/'

# TODO: Take this from impy config, for the respective HImodel 
# Minimal ecm to model hadronic interaction
Ecm_min = 10 # GeV,  10 GeV for SIBYLL23D
Ecm_max = 1e7 # GeV, 1E7 GeV for SIBYLL23D

# Minimal energy for secondaries produced in hadronic interaction
Emin = 10 # GeV

# List of primaries that can undergo hadronic interaction
allowed_primaries = [
    1000010010, # protons only
    # -1000010010, 1000010010, # (anti)protons
    # -1000000010, 1000000010, # (anti)neutrons
]

# List of secondaries produced in hadronic interaction, due to impy not using stable to limit secondaries
allowed_secondaries = [
    22,               # gamma  
    11, -11,          # electron, positron
    211, 111, -211,   # pions
    -2112, 2112,      # (anti)neutron
    -2212, 2212,      # (anti)proton
]
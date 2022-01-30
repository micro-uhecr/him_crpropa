"""Implementation of the hadronic module for CRPropa3

    This implementation is meant as a prototype until the implementation in C++ is finished.

    Date: 9/9/2021
    Leonel Morejon
"""

import sys
from numpy import pi, log, sign, sqrt
from scipy.spatial.transform import Rotation as R
from crpropa import c_light, GeV, Module, ParticleState, Candidate, Vector3d, Random, mass_proton

from config_file import *
sys.path.append(impy_directory)

from impy.kinematics import EventKinematics
from impy import impy_config
from impy.definitions import interaction_model_by_tag, make_generator_instance


# impy_config["user_frame"] = 'laboratory'  # this is not working for some reason


# IMP!!! Currently only a global instance of the generator can be employed1!
# The instance below is used for the moment
# HMRunInstance = make_generator_instance(interaction_model_by_tag['SIBYLL23D'])
# HMRunInstance = make_generator_instance(interaction_model_by_tag['EPOSLHC'])  # 1000 times slower
# HMRunInstance = make_generator_instance(interaction_model_by_tag['PHOJET112'])
# HMRunInstance = make_generator_instance(interaction_model_by_tag['URQMD34'])  # 100 times slower
# HMRunInstance = make_generator_instance(interaction_model_by_tag['PYTHIA8'])  # weird behavior
# HMRunInstance = make_generator_instance(interaction_model_by_tag['QGSJET01C']) 

# mtag = list(interaction_model_by_tag.keys())[19]
HMRunInstance = make_generator_instance(interaction_model_by_tag[mtag])

# Some arbitrary initialization
init_event_kinematics = EventKinematics(
    elab=1e12,    # high enough energy to cover all energies the propagation will use
    p1pdg=2212,
    p2pdg=2212
)
HMRunInstance.init_generator(init_event_kinematics, seed=100001)  # fixed seed for now

class HadronicInteractions(Module):
    '''Prototype class to handle hadronic interactions
    '''
    def __init__(self, matter_density=1e-30, composition={101:100}, distribution=('thermal', 1000), seed=None):
        """The initialization takes as arguments
            - seed           : random numbergenerator seed
            - matter_density : the matter density in units of m-3
            - composition    : a dictionary {pid : particle_density} in arbitrary units
            - distribution   : a list containing (str distribution_type, *distribution_args)
        """
        Module.__init__(self)
        self.matter_density = matter_density  # in m-3
        self.composition = composition
        self.cross_section = 1e-30  # in m2
        self.hadronic_model = HMRunInstance
        self.event_kinematics = init_event_kinematics

        if seed is None:
            self.random_number_generator = Random()  # using th eponymous class from CRPropa
        else:
            self.random_number_generator = Random(seed)
    
    def _compute_interaction_rates(self):
        """Determine the hadronic rates based on inputs: matter density,
        distribution, temperature, and composition.
        """ 
        # ToDo: employ distribution to describe the energy distribution of the target particles
        # ToDo: Compute interaction cross rates based on input parameters

        return self.matter_density * self.cross_section


    def process(self, candidate):
        """This is the function called to operate on candidates (particles)
        and at the moment only works with protons!!!
        """
        # Skip if not a proton!!! or too low energy for the Hadr. Model
        if abs(candidate.current.getId()) not in \
            [ 1000010010, # protons
              1000000010, # neutrons
            ] or candidate.current.getEnergy() < Emin * GeV:
            return
        
        current_step = candidate.getCurrentStep()

        Sigma = self._compute_interaction_rates()

        while current_step > 0:
            # Sampling interaction from the inverse of an exponential distribution
            random_number = self.random_number_generator.rand()
            interaction_step = - log(random_number) / Sigma  # ToDo: optimize sampling for a good range 

            interaction_occurred = current_step > interaction_step
            
        if interaction_occurred:
            candidate.limitNextStep(interaction_step)

            plab = candidate.current.getMomentum().getR()
            g = candidate.current.getLorentzFactor()

            Ecm = mass_proton * sqrt( 2*(g + 1) )
            if Ecm < Ecm_min * GeV:
                return
            else:
                # candidate.limitNextStep(interaction_step)
                plab = candidate.current.getMomentum().getR() / GeV * c_light

            event_kinematics = EventKinematics(
                plab =  plab / GeV * c_light, # projectile momentum, lab frame, GeV/c
                p1pdg = 2212, p2pdg = 2212) # p-p interaction

            secondaries = self.sample_interaction(event_kinematics)
                # Define an arbitrary orthogonal base using the primary direction and arbitrary zenith
                random_angle = 2 * pi * self.random_number_generator.rand()
                vector1, vector2, vector3 = get_orthonormal_base(candidate.current.getDirection(), random_angle)

                for (pid, en, px, py, pz) in secondaries:
                    if pid not in \
                        [     22,             # gamma  
                              11, -11,        # electron, positron
                             211, 111, -211,  # pions
                            -2112, 2112,      # (anti)neutron
                            -2212, 2212,      # (anti)proton
                        ]:
                        continue

                    # Injecting secondaries to CRPropa stack
                    ps = ParticleState()

                    if abs(pid) == 2212:
                        ps.setId(int(sign(pid) * 1000010010)) # crpropa issue with (anti)protons id
                    elif abs(pid) == 2112:
                        ps.setId(int(sign(pid) * 1000000010)) # crpropa issue with (anti)neutrons id
                    else:
                        ps.setId(int(pid))

                    ps.setPosition(candidate.current.getPosition())
                    ps.setEnergy(en * GeV)
                                        
                    # Define orthogonal vector base with arbitrary orientation, but z along the primary direction 
                    vector1 = candidate.current.getDirection().getUnitVector()
                    vector2 = Vector3d(1, 0, 0).cross(candidate.current.getDirection().getUnitVector())
                    vector3 = vector1.cross(vector2)
                    
                    # Rotate the transversal plane a random angle
                    random_angle = 2 * pi * self.random_number_generator.rand()
                    vector2 = vector2.getRotated(vector1, random_angle)
                    vector3 = vector3.getRotated(vector1, random_angle)

                    # Set lengths to secondary momentum components
                    ptot = sqrt(px**2 + py**2 + pz**2)
                    vector1.setR(pz / ptot)
                    vector2.setR(px / ptot)
                    vector3.setR(py / ptot)

                    # Add components to get the resulting momentum of the secondary
                    Secondary_Direction = vector1 + vector2 + vector3
                    Secondary_Direction.setR(1)  # ensure norm
                    
                    ps.setDirection(Secondary_Direction)

                    candidate.addSecondary(Candidate(ps)) # adding secondary to parent's particle stack

                # Remove the covered distance
                current_step -= interaction_step

    def sample_interaction(self, event_kinematics):
        """Calls hadronic model using the interface from impy and 
        the member event_kinematics.
        """

        # TODO: since generating only one event, use method event_generator instead!
        event = list(self.hadronic_model.event_generator(event_kinematics, 1))[0]
        event.filter_final_state()

        # Filtering particles. TODO: implement as input to function
        mask = abs(event.xf) > 0.1

        event_kinematics.boost_cms_to_lab(event)

        secondaries = [sec_properties for sec_properties in 
            zip(event.p_ids[mask], 
                event.en[mask], 
                event.px[mask], 
                event.py[mask], 
                event.pz[mask])]

        return secondaries


def get_orthonormal_base(vector3d, random_angle):
    """Returns a vector orthonormal base where one of the directions
    is aligned to vector3d, and the other two have arbitrary orientation.
    In other words, this function returns three vectors that are normalized, 
    othogonal to each other, and have arbitrary orientation apart from the 
    provided vector3d.
    """

    vector3d.setR(1)
    vector1 = vector3d
    x, y, z = vector1.x, vector1.y, vector1.z
    vector2 = Vector3d(x + y, z + y, 
                       -(x**2 + y**2 + y*(x + z))/z) # orthogonal to vector1
    vector3 = vector1.cross(vector2) # orthogonal to both vector1 and vector2

    # Normalize base
    vector2.setR(1)
    vector3.setR(1)

    # Rotate the transversal plane to vector1 by a random angle
    vector2 = vector2.getRotated(vector1, random_angle)
    vector3 = vector3.getRotated(vector1, random_angle)

        return secondaries
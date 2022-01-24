"""Implementation of the hadronic module for CRPropa3

    This implementation is meant as a prototype until the implementation in C++ is finished.

    Date: 9/9/2021
    Leonel Morejon
"""
import sys
from numpy import log, sign
from crpropa import GeV, Module, ParticleState, Candidate, Vector3d, Random

from config_file import *
sys.path.append(impy_directory)

from impy.kinematics import EventKinematics
from impy import pdata, impy_config
from impy.definitions import interaction_model_by_tag, make_generator_instance


# impy_config["user_frame"] = 'center-of-mass'
impy_config['stable_list'] = []  # no short-lived particles (t < 1h)


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
event_kinematics = EventKinematics(
    elab=1e12,    # high enough energy to cover all energies the propagation will use
    p1pdg=2212,
    p2pdg=2212
)
HMRunInstance.init_generator(event_kinematics)


class HadronicInteractions(Module):
    '''Prototype class to handle hadronic interactions
    '''
    def __init__(self, matter_density=1e-30, composition={101:100}, distribution=('thermal', 1000)):
        """The initialization takes as arguments
            - matter_density : the matter density in units of m-3
            - composition    : a dictionary {pid : particle_density} in arbitrary units
            - distribution   : a list containing (str distribution_type, *distribution_args)
        """
        Module.__init__(self)
        self.matter_density = matter_density
        self.composition = composition
        self.cross_section = 1e-30  # in m2
        self.random_number_generator = Random()  # using th eponymous class from CRPropa
        self.hadronic_model = HMRunInstance
    
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
            [ 1000010010,
              1000000010,
                     111,
                     211
            ] or candidate.current.getEnergy() < 100 * GeV:
            return
        
        current_step = candidate.getCurrentStep()

        Sigma = self._compute_interaction_rates()

        while current_step > 0:
            # Sampling interaction from the inverse of an exponential distribution
            random_number = self.random_number_generator.rand()
            interaction_step = - log(random_number) / Sigma  # ToDo: optimize sampling for a good range 

            interaction_occurred = current_step > interaction_step
            
            if not interaction_occurred:
                return
            else:
                # candidate.limitNextStep(interaction_step)
                # Call to hadronic interaction code and get products
                mass_proton = pdata.mass(2212) # GeV
                secondaries = self.sample_interaction(candidate.current.getEnergy() / GeV + mass_proton)
                # secondaries = self.sample_interaction(12)
                
                for (pid, en, px, py, pz) in secondaries:
                    if pid not in \
                        [     22,        # gamma  
                              11, -11,   # electron, positron
                             211, 111, -211, # pions
                            2112, -2112, # anti-/ neutron
                            2212, -2212, # anti-/ proton
                        ]:
                        continue

                    # Adding secondaries for propagation
                    ps = ParticleState()

                    if abs(pid) == 2212:
                        ps.setId(int(sign(pid) * 1000010010)) # crpropa special treatment of anti-/protons
                    elif abs(pid) == 2112:
                        ps.setId(int(sign(pid) * 1000000010)) # crpropa special treatment of anti-/neutrons
                    else:
                        ps.setId(int(pid))

                    ps.setPosition(candidate.current.getPosition())
                    ps.setEnergy(en * GeV)
                    
                    # Set rotation matrix to bring z into current direction and use to rotate secondaries
                    # Secondary_Direction = candidate.current.getDirection().getRotated(Vector3d(1, 0, 0), 3.14)# ToDo: implement for 3D
                    Secondary_Direction = Vector3d(px, py, pz)
                    Secondary_Direction /= Secondary_Direction.getR()
                    ps.setDirection(Secondary_Direction)

                    candidate.addSecondary(Candidate(ps)) # adding secondary to parent's particle stack

                # Remove the covered distance
                current_step -= interaction_step

    def sample_interaction(self, energy):
        """Calls hadronic model using the interface from impy 

        Arguments:
        ------------
        energy: energy in GeV (lab frame)
        """

        event_kinematics = EventKinematics(
            elab  = energy,
            p1pdg = 2212,
            p2pdg = 2212
        )

        event = list(self.hadronic_model.event_generator(event_kinematics, 1))[0]
        event.filter_final_state()

        # Filtering particles. TODO: implement as input to function     
        mask = abs(event.xf) > 0.1

        secondaries = [sec_properties for sec_properties in 
            zip(event.p_ids[mask], 
                event.en[mask], 
                event.px[mask], 
                event.py[mask], 
                event.pz[mask])]


        return secondaries
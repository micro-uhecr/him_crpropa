"""Implementation of the hadronic module for CRPropa3

    This implementation is meant as a prototype until the implementation in C++ is finished.

    Date: 9/9/2021
    Leonel Morejon
"""

import sys
from numpy import pi, log, sign, sqrt
from scipy.spatial.transform import Rotation as R
from crpropa import c_light, GeV
from crpropa import Module, ParticleState, Candidate, Vector3d, Random

import chromo
from config_file import *
from chromo.models import *
from chromo.kinematics import EventKinematics
from chromo.util import elab2ecm


def get_hi_generator(modelname, initseed):
    """Manages hadronic models, there cannot be 
    two concurrent instances of the same model.
    """
    if modelname not in chromo.models.__all__:
        print('ERROR: Wrong model or model not available.')
        return

    # Some arbitrary initialization
    init_event_kinematics = EventKinematics(
        'proton',
        'proton',
        elab=1e12,    # GeV, high enough energy to cover all energies the propagation will use
    )

    return globals()[modelname](init_event_kinematics, seed=initseed)

HMRunInstance = get_hi_generator('Sibyll23d', 100001)
# allowed_primaries = HMRunInstance.projectiles

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
        # self.cross_section = 1e-30  # in m2

        if seed is None:
            self.random_number_generator = Random()  # using the eponymous class from CRPropa
        else:
            self.random_number_generator = Random(seed)
    
    def _compute_interaction_rates(self, kinematics):
        """Determine the hadronic rates based on inputs: matter density,
        distribution, temperature, and composition.
        """ 
        # ToDo: employ distribution to describe the energy distribution of the target particles
        # ToDo: Compute interaction cross rates based on input parameters

        # sigma = sigma_pp(plab) * 1e-31 # to m2
        sigma = HMRunInstance.cross_section(kinematics).total * 1e-31 # to m2

        # return self.matter_density * self.cross_section
        return self.matter_density * sigma

    def process(self, candidate):
        """This is the function called to operate on candidates (particles)
        and at the moment only works with protons!!!
        """
        projectile_id = candidate.current.getId()

        if abs(projectile_id) == 1000010010:
            projectile_id = int(sign(projectile_id) * 2212)
        elif abs(projectile_id) == 1000000010:
            projectile_id = int(sign(projectile_id) * 2112)

        if projectile_id not in allowed_primaries:
            return
        
        # Set kinematics of of interaction
        # kinematics1 = chromo.kinematics.FixedTarget(1e6 * chromo.constants.TeV, (14, 1), "proton")
        # kinematics2 = chromo.kinematics.CenterOfMass(kinematics1.ecm, 'proton', (14, 1)) # swapping target and projectile
        plab = candidate.current.getMomentum().getR()
        event_kinematics = EventKinematics(
            projectile_id, 'proton', # only p as target for the moment
            plab =  plab / GeV * c_light) # projectile momentum, lab frame, GeV/c

        if not (Ecm_min < event_kinematics.ecm < Ecm_max):
            return

        current_step = candidate.getCurrentStep()
        Sigma = self._compute_interaction_rates(event_kinematics)
        print(f'The interaction rate is {Sigma:3.2e}')

        # Sampling interaction from the inverse of an exponential distribution
        random_number = self.random_number_generator.rand()
        interaction_step = - log(random_number) / Sigma  # ToDo: optimize sampling for a good range 

        interaction_occurred = interaction_step < current_step
        
        if interaction_occurred:
            # Sample colission and get secondaries 
            secondaries = self.sample_interaction(event_kinematics)

            # Arbitrary orthogonal base aligned to primary direction
            random_angle = 2 * pi * self.random_number_generator.rand()
            primary_direction = candidate.current.getDirection()
            vector1, vector2, vector3 = get_orthonormal_base(primary_direction, random_angle)

            # Set interaction position
            step_back = current_step - interaction_step
            xpos = step_back * primary_direction.x
            ypos = step_back * primary_direction.y
            zpos = step_back * primary_direction.z
            interaction_position  = candidate.current.getPosition() - Vector3d(xpos, ypos, zpos)

            Eloss = 0
            for (pid, en, px, py, pz) in secondaries:
                Eloss += en * GeV
                # Injecting secondaries to CRPropa stack
                ps = ParticleState()
                ps.setEnergy(en * GeV)

                if abs(pid) == 2212:
                    ps.setId(int(sign(pid) * 1000010010)) # crpropa issue with (anti)protons id
                elif abs(pid) == 2112:
                    ps.setId(int(sign(pid) * 1000000010)) # crpropa issue with (anti)neutrons id
                else:
                    ps.setId(int(pid))

                ps.setPosition(interaction_position)                
                                    
                # Set lengths to secondary momentum components
                ptot = sqrt(px**2 + py**2 + pz**2)
                cpx = px/ptot * (vector1.x + vector2.x + vector3.x)
                cpy = py/ptot * (vector1.y + vector2.y + vector3.y)
                cpz = pz/ptot * (vector1.z + vector2.z + vector3.z)

                # Add components to get the resulting momentum of the secondary
                Secondary_Direction = Vector3d(cpx, cpy, cpz)
                Secondary_Direction.setR(1)  # enforce unitary norm
                
                ps.setDirection(Secondary_Direction)
                candidate.addSecondary(Candidate(ps)) # adding secondary to parent's particle stack

            # TODO: Can't change primary position (below). why?
            # candidate.current.setPosition(interaction_position)
            candidate.limitNextStep(interaction_step)

            currE = candidate.current.getEnergy()
            candidate.current.setEnergy(currE - Eloss)

    def sample_interaction(self, event_kinematics):
        """Calls hadronic model using the interface from impy and 
        the member event_kinematics.
        """
        HMRunInstance.kinematics = event_kinematics
        
        # TODO: since generating only one event, use method event_generator instead!
        event = list(HMRunInstance(1))[0]

        # TODO: report filter_final_state() does not account for stable config
        event = event.final_state()
        
        # Filtering particles. TODO: implement as input to function
        # mask = (abs(event.xf) > 0.1) * \
        #     (event.en > Emin)

        mask = event.en > 0

        # TODO: Implement substituting not allowed secondaries by its allowed decay products
        secondaries = [sec_properties for sec_properties in 
            zip(event.pid[mask], 
                event.en[mask], 
                event.px[mask], 
                event.py[mask], 
                event.pz[mask]) 
                if sec_properties[0] in allowed_secondaries]

        return secondaries


def get_orthonormal_base(vector3d, random_angle):
    """Returns a vector orthonormal base where one of the directions
    is aligned to vector3d, and the other two have arbitrary orientation.
    The vectors are 
    v1: the vector provided, (x, y, z)
    v2: a vector perpendicular (z-y, x-z, y-x), normalized
    v3: the cross product of v1 and v2

    This method fails in the case of x=y=z=sqrt(3)/3, where the function returns
    v2: sqrt(2)/2, sqrt(2)/2, 0
    v3: the cross product of v1 and v2
    """

    vector1 = vector3d
    vector1.setR(1)
    x, y, z = vector1.x, vector1.y, vector1.z

    norm = sqrt( 2*(1 - x*y - y*z - z*x) )

    if abs(norm) > 1e-6:
        vector2 = Vector3d(z - y, x - z, y - x)  # orthogonal to vector1
        vector2.setR(norm)
    else:
        # only possible when x=y=z=sqrt(3)/3
        vector2 = Vector3d(1/sqrt(2), -1/sqrt(2), 0) # orthogonal to vector1

    vector3 = vector1.cross(vector2) # orthogonal to both vector1 and vector2

    # Rotate the plane transversal to vector1 by a random angle
    vector2 = vector2.getRotated(vector1, random_angle)
    vector3 = vector3.getRotated(vector1, random_angle)

    return vector1, vector2, vector3
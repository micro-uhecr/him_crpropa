"""Implementation of the hadronic module for CRPropa3

    This implementation is meant as a prototype until the implementation in C++ is finished.

    Date: 9/9/2021
    Leonel Morejon
"""

from numpy import pi, log, sign, sqrt, array, cross, arccos, vstack, einsum, vectorize, logical_and
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R
from particle import Particle

from crpropa import Candidate, GeV, Module, ParticleState, Vector3d, Random

from config_file import *
import chromo
from chromo.models import *
from chromo.kinematics import EventKinematics
from chromo.util import elab2ecm, CompositeTarget, EventFrame

# Dealing with proton and neutron pid ambiguity
def pdgid2crpropa(pid):
    """Dealing with provblems due to equivalent
    pid values for protons and neutrons 
    """
    relation = dict([
        (2212, 1000010010), 
        (2112, 1000000010),
        (-2212, -1000010010),
        (-2112, -1000000010) 
    ])
    if pid in relation:
        return relation[pid]
    else:
        return pid

def crpropa2pdgid(pid):
    """Dealing with provblems due to equivalent
    pid values for protons and neutrons 
    """
    relation = dict([
        (1000010010, 2212), 
        (1000000010, 2112),
        (-1000010010, -2212),
        (-1000000010, -2112) 
    ])
    if pid in relation:
        return relation[pid]
    else:
        return pid

def get_hi_generator(modelname, initseed=None):
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

# HMRunInstance = get_hi_generator(mtag, 100001)

class HadronicInteractions(Module):
    '''Prototype class to handle hadronic interactions
    '''
    def __init__(self, matter_density=1e-30, composition={101:100}, distribution=('thermal', 1000), seed=None, Emin=1e6):
        """The initialization takes as arguments
            - seed           : random number generator seed
            - matter_density : the matter density in units of m-3
            - composition    : a dictionary {pid : particle_density} in arbitrary units
            - distribution   : a list containing (str distribution_type, *distribution_args)
        """
        Module.__init__(self)
        self.matter_density = matter_density  # in m-3
        self.composition = composition
        self.Emin = Emin

        if seed is None:
            self.random_number_generator = Random()  # using the eponymous class from CRPropa
        else:
            self.random_number_generator = Random(seed)

        self.hi_engine = get_hi_generator(mtag, seed)
    
    def _compute_interaction_rates(self, kinematics):
        """Determine the hadronic rates based on inputs: matter density,
        distribution, temperature, and composition.
        """ 
        # ToDo: employ distribution to describe the energy distribution of the target particles
        # ToDo: Compute interaction cross rates based on input parameters

        # sigma = sigma_pp(plab) * 1e-31 # to m2
        sigma = self.hi_engine.cross_section(kinematics).total * 1e-31 # to m2

        # return self.matter_density * self.cross_section
        return self.matter_density * sigma

    def process(self, candidate):
        """This is the function called to operate on candidates (particles)
        and at the moment only works with protons!!!
        """
        target = crpropa2pdgid(1000010010)
        projectile = crpropa2pdgid(candidate.current.getId())
        currE = candidate.current.getEnergy() / GeV # in GeV
        
        m1 = Particle.from_pdgid(projectile).mass / 1e3 # in GeV
        m2 = Particle.from_pdgid(target).mass / 1e3 # in GeV
        Ecm = elab2ecm(currE, m1, m2)
        evkin = chromo.kinematics.CenterOfMass(Ecm, target, projectile)

        # Invert particle order if projectile not allowed by generator
        if projectile not in self.hi_engine.projectiles:
            if (target in self.hi_engine.projectiles) and (projectile in self.hi_engine.targets):
                # swapping target and projectile with equivalent c.m. energy
                evkin = chromo.kinematics.CenterOfMass(Ecm, target, projectile)
            else:
                # Interaction not possible with chosen generator
                return

        # Set kinematics of of interaction
        try:
            self.hi_engine.kinematics = evkin
        except:
            # Avoiding all exceptions from generator due to kinematics (Emin, pid unavailable, etc.)
            # Skipping step sampling if kinematics are not supported
            return

        current_step = candidate.getCurrentStep()
        Sigma = self._compute_interaction_rates(evkin)

        # Sampling interaction from the inverse of an exponential distribution
        random_number = self.random_number_generator.rand()
        interaction_step = - log(random_number) / Sigma  # ToDo: optimize sampling for a good range 

        interaction_occurred = interaction_step < current_step
        
        if interaction_occurred:
            # Sample colission and get secondaries 
            pids, energies, momenta = self.sample_interaction()
            Eloss = energies.sum() * GeV
            primary_direction = array(candidate.current.getDirection())

            # Set interaction position in CRPropa corrdinate system
            step_back = (current_step - interaction_step) * primary_direction
            interaction_position  = candidate.current.getPosition() - Vector3d(step_back[0], step_back[1], step_back[2])

            # Arbitrary orthogonal base aligned to primary direction
            random_phi = 2 * pi * self.random_number_generator.rand()            
            arbitrary_phi_rotation = R.from_euler('z', random_phi).as_matrix()
            rotation_axis = cross(array([0, 0, 1]), primary_direction)
            theta = arccos(primary_direction[2]) # in radians
            z_alignment_to_primary_direction = R.from_rotvec(theta * rotation_axis).as_matrix()
            transformation = arbitrary_phi_rotation.dot(z_alignment_to_primary_direction)

            momenta = momenta.dot(transformation)

            for pid, en, pvector in zip(pids, energies, momenta):
                # Injecting secondaries, adding secondary to parent's particle stack
                direction = Vector3d(pvector[0], pvector[1], pvector[2])
                ps = ParticleState(int(pid), en * GeV, interaction_position, direction)
                candidate.addSecondary(Candidate(ps))

            # TODO: improve limiting step condition based on cross section
            candidate.limitNextStep(interaction_step)

            # TODO: Verify if primary should be deactivated
            # candidate.current.setEnergy(currE - Eloss)
            candidate.current.setEnergy(0)

    def sample_interaction(self):
        """Calls hadronic model and returns products
           Returns:
           - particles ids
           - energies [in GeV]
           - momenta as unitary vectors 
        """        
        # TODO: since generating only one event, use method event_generator instead!
        event = list(self.hi_engine(1))[0]

        event = event.final_state()
        
        # Applying boost to lab frame
        event.kin.apply_boost(event, EventFrame.FIXED_TARGET)

        # Filtering particles
        mask = logical_and(event.en > self.Emin, [pid in allowed_secondaries for pid in event.pid])
        
        energies = event.en[mask]
        momenta = vstack([event.px[mask], event.py[mask], event.pz[mask]]).T
        momenta = einsum('ij, i -> ij', momenta, 1 / norm(momenta, axis=1)) # normalizing
        pids = [pdgid2crpropa(pid) for pid in event.pid[mask]] # Fixing some pid values

        # TODO: Implement substituting not allowed secondaries by its allowed decay products

        # return (pids[:4], event.en[:4], momenta[:4, :])
        return (pids, energies, momenta)
        # return (array([1000260560, 1000260560]), array([.1, .1]), array([[0, 0, 1], [0, 0, 1]]))

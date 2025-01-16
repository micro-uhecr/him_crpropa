"""Implementation of the hadronic module for CRPropa3

    This implementation is meant as a prototype until the implementation in C++ is finished.

    Date: 9/9/2021
    Leonel Morejon
    
    Updated by Julien Dörner on 16/01/2025
"""

__version__ = "0.0.1"

from numpy import pi, log, sqrt, array, cross, arccos, vstack, einsum, logical_and, logspace
import numpy as np
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import interp1d
from particle import Particle

from crpropa import Candidate, mass_proton, c_light, GeV, Module, ParticleState, Vector3d, Random, ConstantDensity

from him_crpropa.config_file import *
import chromo
from chromo.models import *
from chromo.kinematics import CenterOfMass
from chromo.kinematics import EventKinematicsWithRestframe as EventKinematics
from chromo.util import elab2ecm, CompositeTarget, EventFrame

mp = mass_proton * c_light**2 / GeV

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


def sample_model_xsec(hi_generator):
    """Sampling the model cross section for proton-proton 
    interactions, as a function of the laboratory momentum 
    plab in GeV. Returns xsec in milibarn.
    """
    csec_vals = []
    ecm_vals = logspace(1, 6, 30)
    plab_vals = []
    for ecm in ecm_vals:
        kin = CenterOfMass(ecm, "p", "p")
        plab_vals.append(sqrt(kin.elab**2 - (kin.p1.A * mp)**2))
        csec_vals.append(hi_generator.cross_section(kin).inelastic)
        
    return interp1d(plab_vals, csec_vals)

def sigma_pp(plab):
    """Cross section for proton-proton interactions based on the PDG fit, as
    a function of the laboratory momentum plab in GeV. Returns xsec in milibarn.
    The inelastic cross section is modeled as a fraction of this total cross section.

    Reference: C. Patrignani 2016 Chinese Phys. C 40 100001
    """
    mp = 0.938272 # GeV
    M = 2.1206 # GeV
    H = 0.272 # mb
    P, R1, R2 = 34.41, 13.07, 7.394 # in mb
    eta1, eta2 = 0.4473, 0.5486 # dimenssionless

    ecm2 = 2*(mp**2 + mp*sqrt(plab**2 + mp**2)) # GeV
    sab = (2*mp + M)**2 # GeV
    xsec_total = H * log(ecm2/sab)**2 + P + R1*(ecm2/sab)**-eta1 - R2*(ecm2/sab)**-eta2

    return 0.81 * xsec_total

class HadronicInteractions(Module):
    '''Prototype class to handle hadronic interactions
    '''
    def __init__(self, matter_density=1e-30, composition={101:100}, distribution=('thermal', 1000), seed=None, Emin=1e6, mtag='Sibyll23d', model_xsec=False):
        """The initialization takes as arguments
            - seed           : random number generator seed
            - matter_density : the matter density in units of m-3 (default) or an instance of Density defined in CRPropa
            - composition    : a dictionary {pid : particle_density} in arbitrary units
            - distribution   : a list containing (str distribution_type, *distribution_args)
        """
        Module.__init__(self)

        self.matter_density = matter_density
        self.composition = composition
        self.Emin = Emin
        self.allowed_secondaries = allowed_secondaries
        
        self.model_xsec = model_xsec
        self._hi_engine = object
        self.hi_engine = mtag

        if seed is None:
            self.random_number_generator = Random()  # using the eponymous class from CRPropa
        else:
            self.random_number_generator = Random(seed)
        
        self.limit_secondaries()

    @property
    def hi_engine(self):
        return self._hi_engine

    @hi_engine.setter
    def hi_engine(self, new_engine):
        """Setting or changing the hadronic interactions generator.
        Allowed values are any of the mtag values available in chromo
        or directly the result of get_hi_generator.
        """

        if type(new_engine) == str:
            if new_engine in str(type(self._hi_engine)):
                return
            self._hi_engine = get_hi_generator(new_engine)
        else:
            if new_engine == self._hi_engine:
                return
            self._hi_engine = new_engine
        
        if self.model_xsec:
            print('Sampling interaction cross sections.')
            self.xsec = sample_model_xsec(self._hi_engine)
            print('Sampling interaction cross sections completed.')
        else:
            self.xsec = sigma_pp

        self.__name__ = f'HIM ({type(new_engine)})'
        
    @property
    def matter_density(self):
        return self._matter_density
    
    @matter_density.setter
    def matter_density(self, matter_density):
        """Setting the matter density: if a number, constant density in units of m^-3,
        otherwise takes an instance of density class from CRPropa.
        """
        if type(matter_density) in [int, float, np.int64, np.float64]:
            self._matter_density = ConstantDensity(matter_density, 0, 0) # in m-3
        else:
            self._matter_density = matter_density # instance based on class Density

    def limit_secondaries(self, max_lifetime=1e30):
        """Restricts the secondaries to those with a 
           lifetime greater than max_lifetime in seconds.
        """
    
        for part in Particle.all():            
            if (part.lifetime is not None) and (part.lifetime < max_lifetime * 1e9):
                try:
                    self.hi_engine.set_unstable(part.pdgid)
                    print(part.name)
                except:
                    print(f"WARNING: {part.name} not allowed as secondary.")
                
    
    def compute_interaction_rates(self, kinematics, matter_density):
        """Determine the hadronic rates based on inputs: matter density,
        distribution, temperature, and composition.
        """ 
        # ToDo: employ distribution to describe the energy distribution of the target particles
        # ToDo: Compute interaction cross rates based on input parameters

        sigma = self.xsec(kinematics.plab) * 1e-31 # to m2
        # sigma = self.hi_engine.cross_section(kinematics).total * 1e-31 # to m2

        return matter_density * sigma

    def process(self, candidate):
        """This is the function called to operate on candidates (particles)
        and at the moment only works with protons!!!
        """
        if not candidate.current.getId() == 1000010010:
            return
        
        target = crpropa2pdgid(1000010010)
        projectile = crpropa2pdgid(candidate.current.getId())
        currE = candidate.current.getEnergy() / GeV # in GeV
        matter_density = self._matter_density.getNucleonDensity(candidate.current.getPosition())
        
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
        Sigma = self.compute_interaction_rates(evkin, matter_density)
        if Sigma == 0:
            return

        # Sampling interaction from the inverse of an exponential distribution
        random_number = self.random_number_generator.rand()
        interaction_step = - log(random_number) / Sigma  # ToDo: optimize sampling for a good range 

        interaction_occurred = interaction_step < current_step
        
        if interaction_occurred:
            # Sample colission and get secondaries 
            pids, energies, momenta = self.sample_interaction()
            Eloss = energies.sum() * GeV
            crpropa_direction = candidate.current.getDirection()
            primary_direction = array([crpropa_direction.x, crpropa_direction.y, crpropa_direction.z])
        
            
            # Set interaction position in CRPropa corrdinate system
            interaction_position = candidate.current.getPosition() - candidate.current.getDirection() * interaction_step

            # Arbitrary orthogonal base aligned to primary direction
            random_phi = 2 * pi * self.random_number_generator.rand()            
            arbitrary_phi_rotation = R.from_euler('z', random_phi).as_matrix()
            primary_direction_xy = primary_direction.copy()
            primary_direction_xy[2] = 0
            rotation_axis = cross(array([0, 0, 1]), primary_direction_xy)
            theta = arccos(primary_direction[2]) # in radians
            z_alignment_to_primary_direction = R.from_rotvec(theta * rotation_axis).as_matrix()
            transformation = arbitrary_phi_rotation.dot(z_alignment_to_primary_direction)

            momenta = momenta.dot(transformation)

            for pid, en, pvector in zip(pids, energies, momenta):
                # Injecting secondaries, adding secondary to parent's particle stack
                direction = Vector3d(pvector[0], pvector[1], pvector[2])
                ps = ParticleState(int(pid), en * GeV, interaction_position, direction)
                candidate.addSecondary(Candidate(ps))

            # TODO: improve limiting step condition should be done if interaction did not occur
            candidate.limitNextStep(interaction_step)

            # TODO: Primary should be deactivated
            # candidate.current.setEnergy(currE - Eloss)
            candidate.setActive(False)

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
        mask = event.en > self.Emin
        mask = logical_and(event.en > 0, [pid in self.allowed_secondaries for pid in event.pid])
        
        energies = event.en[mask]
        momenta = vstack([event.px[mask], event.py[mask], event.pz[mask]]).T
        momenta = einsum('ij, i -> ij', momenta, 1 / norm(momenta, axis=1)) # normalizing
        pids = [pdgid2crpropa(pid) for pid in event.pid[mask]] # Fixing some pid values

        # TODO: Implement substituting not allowed secondaries by its allowed decay products

        # return (pids[:4], event.en[:4], momenta[:4, :])
        return (pids, energies, momenta)
        # return (array([1000260560, 1000260560]), array([.1, .1]), array([[0, 0, 1], [0, 0, 1]]))

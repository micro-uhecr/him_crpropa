from crpropa import *
from him_crpropa.hadronic_module import * 
import sys

def TestRun1D(N=100, filename='Test1.txt', density=1e13, seed=None):
    """Testing function for hadronic module in 1D

    Primary protons of 1 EeV, one directionalin 1 pc spherical volume.
    Cutoff energy is 1 TeV, only hadronic interactions for protons with
    Cross section 10mb

    """

    #specifying the output
    particle_positions = TextOutput(filename, Output.Event1D)
    
    # Detection grid

    sim = ModuleList()
    det = Observer()
    det.add(ObserverSurface(Sphere(Vector3d(0,0,0), 1 * pc)))
    det.onDetection(particle_positions)
    sim.add(det)

    module1 = SimplePropagation(0.01*pc, 0.5*pc)
    module3 = HadronicInteractions(matter_density=density, seed=seed)
    sim.add(module1)
    sim.add(module3)
    sim.add(SphericalBoundary(Vector3d(0, 0, 0), 1.01 * pc))
    
    sim.add( MinimumEnergy( 100 * TeV) )    
    
    # Specifying the source
    source = Source()
    source.add(SourceParticleType(nucleusId(1, 1)))
    source.add(SourceDirection(Vector3d(0, 0, 1)))
    source.add(SourcePosition(0))
    source.add(SourceEnergy(1 * EeV))

    # Run particle
    sim.setShowProgress(True)
    sim.run(source, N)
    
def TestRun1D_tracks(N=100, filename='Test1.txt', density=1e13, seed=None):
    """Testing function for hadronic module in 1D

    Primary protons of 1 EeV, one directionalin 1 pc spherical volume.
    Cutoff energy is 1 TeV, only hadronic interactions for protons with
    Cross section 10mb

    """

    #specifying the output
    particle_positions = TextOutput(filename, Output.Trajectory3D)
    
    # Detection grid

    sim = ModuleList()
    det = Observer()
    # det.add(ObserverTracking(Vector3d(0,0,0), pc, 0.01*pc))
    det.add(ObserverSurface(Sphere(Vector3d(0,0,0), 1 * pc)))
    det.onDetection(particle_positions)
    sim.add(det)
        
    module1 = SimplePropagation(0.3*pc, 1.5*pc)
    if seed == None:
        module3 = HadronicInteractions(matter_density=density)
    else:
        module3 = HadronicInteractions(matter_density=density, seed=seed)

    sim.add(module1)
    sim.add(module3)
    sim.add(SphericalBoundary(Vector3d(0, 0, 0), 1.01 * pc))
    
    sim.add( MinimumEnergy( 10 * TeV) )    
    
    # Specifying the source
    source = Source()
    source.add(SourceParticleType(nucleusId(1, 1)))
    source.add(SourceDirection(Vector3d(0, 0, 1)))
    source.add(SourcePosition(0))
    source.add(SourceEnergy(1 * EeV))

    # Run particle
    sim.setShowProgress(True)
    sim.run(source, N)

def TestReproducibility(N=100, filename='Test1.txt', density=1e13):
    """Testing that the same seeds yield unchanged results
    """
    # Fixing CRPropa seed
    Random_seedThreads(1987)
    # Fixing HM seed
    
    
    #specifying the output
    particle_positions = TextOutput(filename, Output.Trajectory3D)
    
    # Detection grid

    sim = ModuleList()
    det = Observer()
    
    det.add(ObserverSurface(Sphere(Vector3d(0,0,0), 1 * pc)))
    det.onDetection(particle_positions)
    sim.add(det)
        
    module1 = SimplePropagation(0.3*pc, 1.5*pc)
    module3 = HadronicInteractions(matter_density=density, seed=1298)

    sim.add(module1)
    sim.add(module3)
    sim.add(SphericalBoundary(Vector3d(0, 0, 0), 1.01 * pc))
    
    sim.add( MinimumEnergy( 1 * PeV) )    
    
    # Specifying the source
    source = Source()
    source.add(SourceParticleType(nucleusId(1, 1)))
    source.add(SourceDirection(Vector3d(0, 0, 1)))
    source.add(SourcePosition(0))
    source.add(SourceEnergy(50 * TeV))

    # Run particle
    sim.setShowProgress(True)
    sim.run(source, N)

def TestOneInteractionLength(N=100, filename='Test1.txt', energy=1, seed=None):
    """Testing function for hadronic module isotropic source.

    Primary protons of given energy (EeV), isotropic source, in spherical volume
    of 1 pc and density such that the radius equals one interaction length.
    No cutoff energy, only hadronic interactions for protons with cross section
    given by PDG model.
    """

    def sigma_pp(plab):
            """Cross section for proton-proton interactions based on the PDG fit, as
            a function of the laboratory momentum plab in GeV.

            Reference: C. Patrignani 2016 Chinese Phys. C 40 100001
            """
            mp = 0.938272 # GeV
            M = 2.1206 # GeV
            H = 0.272 # mb
            P, R1, R2 = 34.41, 13.07, 7.394 # in mb
            eta1, eta2 = 0.4473, 0.5486 # dimenssionless

            ecm2 = 2*(mp**2 + mp*sqrt(plab**2 + mp**2)) # GeV
            sab = (2*mp + M)**2 # GeV

            xsec = H * log(ecm2/sab)**2 + P + R1*(ecm2/sab)**-eta1 - R2*(ecm2/sab)**-eta2

            return xsec 

    density =  1 / sigma_pp(energy * 1e9) / 1e-31 / pc # m^-3

    #specifying the output
    particle_positions = TextOutput(filename, Output.Event1D)
    
    # Detection grid

    sim = ModuleList()
    det = Observer()
    det.add(ObserverSurface(Sphere(Vector3d(0,0,0), 1 * pc)))
    det.onDetection(particle_positions)
    sim.add(det)

    module1 = SimplePropagation(0.01*pc, 0.5*pc)
    module3 = HadronicInteractions(matter_density=density, seed=seed)
    sim.add(module1)
    sim.add(module3)
    sim.add(SphericalBoundary(Vector3d(0, 0, 0), 1.01 * pc))   
    
    # Specifying the source
    source = Source()
    source.add(SourceParticleType(nucleusId(1, 1)))
    source.add(SourceIsotropicEmission())
    source.add(SourceEnergy(1 * EeV))

    # Run particle
    sim.setShowProgress(True)
    sim.run(source, N)


if __name__ == "__main__":
    # TestRun1D(100000, 'output.txt')
    Nprimaries = 100
    output_filename = 'output.txt'
    if len(sys.argv) > 2:
        Nprimaries = int(sys.argv[1])
    if len(sys.argv) == 3:
        output_filename = sys.argv[2]

    Random_seedThreads(1987)
    # TestRun1D(Nprimaries, output_filename, seed=1967)
    # TestRun1D_tracks(Nprimaries, output_filename, seed=1967, density=1e22)
    # TestReproducibility(Nprimaries, output_filename, density=1e19)

    # print(Random_seedThreads.getSeed())
    TestOneInteractionLength(Nprimaries, output_filename)

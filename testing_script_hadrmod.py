from crpropa import *
from hadronic_module import * 

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

if __name__ == "__main__":
    # TestRun1D(100000, 'output.txt')
    Nprimaries = int(sys.argv[1])
    output_filename = sys.argv[2]

    Random_seedThreads(1987)
    # TestRun1D(Nprimaries, output_filename, seed=1967)
    TestRun1D_tracks(Nprimaries, output_filename, seed=1967, density=1e22)
    # TestReproducibility(Nprimaries, output_filename, density=1e19)

    # print(Random_seedThreads.getSeed())

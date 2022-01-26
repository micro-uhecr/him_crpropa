from crpropa import *

import sys
sys.path.append('CRPropa_works')  # location of hadronic module

from hadronic_module import * 

def TestRun1D(N=100, filename='Test1.txt', density=1e13):
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
    module3 = HadronicInteractions(matter_density=density)
    sim.add(module1)
    sim.add(module3)
    sim.add(SphericalBoundary(Vector3d(0, 0, 0), 1.01 * pc))
    
    sim.add( MinimumEnergy( 1 * TeV) )    
    
    # Specifying the source
    source = Source()
    source.add(SourceParticleType(nucleusId(1, 1)))
    source.add(SourceDirection(Vector3d(0, 0, 1)))
    source.add(SourcePosition(0))
    source.add(SourceEnergy(1 * EeV))

    # Run particle
    sim.setShowProgress(True)
    sim.run(source, N)
    
    
if __name__ == "__main__":
    # TestRun1D(100000, 'output.txt')
    Nprimaries = int(sys.argv[1])
    output_filename = sys.argv[2]

    TestRun1D(Nprimaries, output_filename)

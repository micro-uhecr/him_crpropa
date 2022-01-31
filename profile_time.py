import numpy as np
from tqdm import tqdm
import time
from testing_script_hadrmod import *

hmodel_list = ["PHOJET112", "SIBYLL23D", "EPOSLHC", "QGSJET01C", "URQMD34"]

def base_profile(Nsim=10, Np=100, density=1e17):
    """Runs all models a number of cycles and reports on the 
    average and std for time duration.
    Nsim : number of cycles
    Np : number of primaries per cycle
    """

    tprofile = np.zeros((Nsim, len(hmodel_list)))

    with tqdm(total=100) as pbar:
        for nhm, hmodel in enumerate(hmodel_list):

            HMRunInstance = get_hi_generator(hmodel)

            for simnum in range(1, Nsim):

                Nprimaries = 10
                output_filename = f"output_{hmodel}_Np{Nprimaries:}_{simnum:02d}.txt"

                Random_seedThreads(1987)

                tstart = time.time()
                TestRun1D_tracks(Nprimaries, output_filename, seed=None, density=density)
                dt = time.time() - tstart

                tprofile[simnum, nhm] = dt

                pbar.update(100 * Nsim/len(hmodel_list))

    np.savetxt('time_profile_' + ('_'.join(hmodel_list) + '.txt'), tprofile)


if __name__ == "__main__":

    base_profile()


    

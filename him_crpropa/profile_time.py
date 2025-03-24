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
                output_filename = f"output_{hmodel}_Np{Np:}_{simnum:02d}.txt"

                tstart = time.time()
                TestRun1D_tracks(Np, output_filename, seed=None, density=density)
                dt = time.time() - tstart

                tprofile[simnum, nhm] = dt

                pbar.update(100 * Nsim/len(hmodel_list))

        np.savetxt('time_profile_' + ('_'.join(hmodel_list) + '.txt'), tprofile)


if __name__ == "__main__":
    import os

    density_grid = np.logspace(15, 30, 30, endpoint=False)

    for d in density_grid[:]:
        dirname = f'density={d:3.2e}'
        os.makedirs(dirname, exist_ok=True)

        base_profile(Nsim=10, Np=1000)

        for filename in os.listdir('./'):
            if filename.startswith('output') and filename.endswith('txt'):
                os.rename(filename, os.path.join(dirname, filename))

        tprof_filename = 'time_profile_' + ('_'.join(hmodel_list) + '.txt')
        os.rename(tprof_filename, os.path.join(dirname, tprof_filename))


    

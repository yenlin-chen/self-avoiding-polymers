import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from os import makedirs

from polymer import RandomWalkPolymer_xy

if __name__ == '__main__':

    data_dir = 'data'

    makedirs(data_dir)

    ####################################################################
    # define parameters
    ####################################################################

    # this will generate <n_iteration> .npy files
    # each with <n_samples> chains
    n_samples = 1000
    n_iterations = 5

    # parameters
    target_length_ls = [25, 50, 75, 100, 125, 150, 175, 200]

    ####################################################################
    # do sampling
    ####################################################################

    for iter_idx in range(n_iterations):
        print(f'Iteration {iter_idx+1} of {n_iterations}')
        for target_length in target_length_ls:

            # instantiate random walker
            walker = RandomWalkPolymer_xy()

            # define placeholders for sampling results
            polymers = np.empty((n_samples, target_length, 2))
            coms = np.empty((n_samples, 2)) # center of mass
            polymers_com = np.empty(polymers.shape) # translate by com

            # sample from population
            for sample_idx in tqdm(range(n_samples),
                                   desc=(f' >> polymer length '
                                         f'{target_length}'),
                                   ascii=True, dynamic_ncols=True):
                polymer = walker.grow_to(target_length)
                polymers[sample_idx] = polymer

            np.save(f'{data_dir}/len_{target_length}-{n_samples}_'
                    f'samples-{iter_idx}',
                    polymers)

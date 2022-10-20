import numpy as np
import matplotlib.pyplot as plt
from os import makedirs, listdir
import imageio

if __name__ == '__main__':

    data_dir = 'data'
    plot_dir = 'figures'
    cache_dir = 'growth'

    makedirs(f'{plot_dir}/{cache_dir}', exist_ok=True)

    ####################################################################
    # define parameters
    ####################################################################

    # number of samples to include in .gif
    n_samples = 20

    # polymer length to plot
    target_length = 100

    ####################################################################
    # read data
    ####################################################################

    file = [f for f in listdir(data_dir) if f'len_{target_length}' in f][0]

    polymers = np.load(f'{data_dir}/{file}')

    coms = np.ma.average(polymers, axis=1)
    polymers_com = polymers - coms[:,np.newaxis,:]

    ####################################################################
    # plot and save all frames
    ####################################################################

    # parameters for plotter
    xylim = int(target_length/3)+5

    for n_idx in range(n_samples):
        for idx in range(target_length):

            fig, ax1 = plt.subplots(1,1, figsize=(8,8))

            ax1.axis('equal')
            ax1.set_axisbelow(True)
            ax1.set_aspect('equal', 'box')
            ax1.grid()

            ax1.set_xlim(-xylim,xylim)
            ax1.set_ylim(-xylim,xylim)
            ax1.tick_params(which='both', length=0)

            # plot grown chains
            if n_idx > 0:
                ax1.plot(*polymers[:n_idx].T, alpha=0.5, c='C0', zorder=7)
            # plot growing chains
            ax1.plot(*polymers[n_idx,:idx+1].T,
                     alpha=1, c='m', zorder=8)
            # plot head of chain
            ax1.scatter(*polymers[n_idx,idx].T,
                        alpha=1, c='C1', zorder=10)

            # mark origin
            ax1.scatter(0, 0, c='cyan', zorder=9, marker='s',
                        label='head of chain')
            # mark end of grown chains
            unique_endings, counts_endings = np.unique(
                polymers[:n_idx,-1,:],
                axis=0,
                return_counts=True
            )
            for c_idx in np.unique(counts_endings):
                ax1.scatter(*unique_endings[counts_endings==c_idx].T,
                           c=f'C{c_idx+2}', zorder=9, alpha=1,
                           label=f'tail (duplicity: {c_idx})')
            ax1.legend(loc='upper right')

            # save images to disk
            plt.savefig(f'{plot_dir}/{cache_dir}/samples-len{target_length}-'
                        f'{n_idx*target_length+idx}.png',
                        dpi=150, bbox_inches='tight')
            plt.close()

    ####################################################################
    # generate gif
    ####################################################################

    n_frames = n_idx*target_length+idx+1

    images = []

    for idx in range(n_frames):
        file_path = f'{plot_dir}/{cache_dir}/samples-len{target_length}-{idx}.png'
        images.append(imageio.v2.imread(file_path))

    # Make it pause at the end so that the viewers can ponder
    for _ in range(10):
        images.append(imageio.v2.imread(file_path))

    imageio.mimsave(f'{plot_dir}/growth-len{target_length}.gif', images,
                    duration=0.01)

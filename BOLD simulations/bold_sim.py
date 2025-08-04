import numpy as np
import matplotlib.pyplot as plt
from neurolib.models.aln import ALNModel
from neurolib.utils.loadData import Dataset
from scipy.stats import pearsonr
import pandas as pd
import sys
import os
import neurolib
print(neurolib.__file__)

run_example = False
run_series = False
run_matrix = False
save_conn = True

# set simulation parameters
pre_duration = 5000
sim_duration = 60000

# check that python is running in neurolib environment
print(f"Python executable: {sys.executable}")  # Should point to neurolib_env
print(f"Running in: {sys.prefix}")  # Should show neurolib_env path

#set random seed
random_seed = 1
np.random.seed(random_seed)

# Create output directory if it doesn't exist
output_dir = "simulation_outputs"
os.makedirs(output_dir, exist_ok=True)

# load a dataset
ds = Dataset("gw")

# initialize the model with structural connectivity and delay matrices
model = ALNModel(Cmat=ds.Cmat, Dmat=ds.Dmat, seed=random_seed)

# set parameters
params = {
    'duration': sim_duration,      # simulation duration, in ms
    'bold_dt': 100,         # bold sampling time step, in ms (isn't it actually 2000 ms?)
    'dt': 0.1,              # simulation time step, in ms
    'bold': True,
    'mue_ext_mean': 1.6,    #
    'mui_ext_mean': 1.6,    #
    'sigma_ou': 0.1,        #
    'seed': random_seed,              # Main seed
    'sampling_seed': random_seed,     # BOLD subsampling seed
    'noise_seed': random_seed,        # Noise generation seed
    'chaos_seed': random_seed,        # Initial condition seed
    'random_seed': random_seed,       # Fallback for other RNG
    #'Ke_gl': 10.0,         # global coupling (required for multi-region)
    'bold_rho': 0.34,
    'bold_alpha': 0.32,
    'bold_V0': 0.02,
    'bold_k1': 3.72,
    'bold_k2': 0.53,
    'bold_k3': 0.53,
    'bold_Gamma': 0.41,
    'bold_K': 0.65,
    'bold_Tau': 0.98,
}
model.params.update(params)

if save_conn:
    #save connectivity matrix
    np.savetxt('gw_connectivity_matrix.csv', ds.Cmat, delimiter=',')

if run_example:
    #run simulation
    model.run()

    # Downsample neural activity to 1 ms time steps
    downsample_factor = 10  # Since original dt=0.1 ms, 10x downsampling gives 1 ms
    downsampled_neural = model.rates_exc[:, ::downsample_factor]  # Take every 10th sample
    downsampled_t = model.t[::downsample_factor]  # Corresponding downsampled time points

    # Save downsampled neural activity as CSV
    neural_df = pd.DataFrame(downsampled_neural.T,
                            columns=[f"Region_{i}" for i in range(downsampled_neural.shape[0])],
                            index=downsampled_t)
    neural_df.index.name = "Time_ms"
    neural_df.to_csv(os.path.join(output_dir, f"neural_activity_seed_{random_seed}.csv"))

    # Save BOLD signals for each node as CSV
    bold_signals = model.BOLD.BOLD.T  # Transpose to (timepoints, regions)
    bold_df = pd.DataFrame(bold_signals,
                          columns=[f"Region_{i}" for i in range(bold_signals.shape[1])],
                          index=model.BOLD.t_BOLD)
    bold_df.index.name = "Time_s"
    bold_df.to_csv(os.path.join(output_dir, f"bold_signals_seed_{random_seed}.csv"))

    # verify
    print(f"Downsampled neural activity shape: {downsampled_neural.shape}")
    print(f"BOLD data shape: {bold_signals.shape}")
    print(f"Simulation duration: {model.t[-1]} ms")

    # plot BOLD for 1st 5 regions
    print("Plotting BOLD")
    plt.figure(figsize=(12, 6))
    for region in range(5):  # First 5 regions
        plt.plot(model.BOLD.t_BOLD, model.BOLD.BOLD[region, :],
                 'o-',
                 label=f'Region {region+1}')
    plt.xlabel('Time (s)')
    plt.xlim(pre_duration, sim_duration) #exclude pre duration
    plt.title(f'BOLD Signals (2000ms intervals)')
    plt.legend(bbox_to_anchor=(1.05, 1))  # Move legend outside
    plt.tight_layout()
    plt.show()

    # plot excitatory firing rates for 1st 5 regions
    print("Plotting excitatory firing rates")
    plt.figure(figsize=(12, 6))
    for i in range(5):  # num regions
        plt.plot(downsampled_t, downsampled_neural[i], label=f'Region {i+1} (Exc)')
    plt.xlabel('Time (ms)')
    plt.xlim(pre_duration, sim_duration)
    plt.ylabel('Firing rate (Hz)')
    plt.title('Excitatory Population Firing Rates')
    plt.legend()
    plt.show()


    # Compute FC matrix (Pearson correlations)
    print("Plotting BOLD FC matrix")
    bold_signals = model.BOLD.BOLD.T  # Shape: (timepoints, regions)
    n_regions = bold_signals.shape[1]
    fc_matrix = np.zeros((n_regions, n_regions))
    for i in range(n_regions):
        for j in range(n_regions):
            fc_matrix[i,j] = pearsonr(bold_signals[:,i], bold_signals[:,j])[0]
    # Plot FC matrix
    plt.figure(figsize=(8, 6))
    plt.imshow(fc_matrix, cmap='coolwarm', vmin=0, vmax=1)
    plt.colorbar(label='Functional Connectivity (r)')
    plt.title('BOLD Functional Connectivity Matrix')
    plt.xlabel('Region')
    plt.ylabel('Region')
    plt.show()

    global_fc = np.mean(np.abs(fc_matrix[np.triu_indices_from(fc_matrix, k=1)]))
    print(f"Global functional connectivity strength: {global_fc:.4f}")



if run_series:
    #run many simulations and calculate fc matrix
    num_simulations = 1
    results = []
    n_regions = 80
    header = ["sim", "global_fc"] + [f"fc_{i}" for i in range(n_regions * n_regions)]
    with open('fc_results.csv', 'w') as f:
        # Write header
        f.write(",".join(header) + "\n")

        for seed in range(num_simulations):
            # Run simulation and compute FC matrix
            model.params.update({'seed': seed})
            model.randomState = np.random.RandomState(seed)
            if hasattr(model, 'noise_input'):
               model.noise_input = []  # Clear noise buffers
            model.run()
            fc_matrix = np.corrcoef(model.BOLD.BOLD.T, rowvar=False)
            global_fc = np.mean(np.abs(fc_matrix[np.triu_indices_from(fc_matrix, k=1)]))
            print(f"sim {seed} global fc: {global_fc:.4f}")

            # Prepare row: sim#, global_fc, then flattened FC values
            row = [str(seed), f"{global_fc:.4f}"] + [f"{x:.4f}" for x in fc_matrix.flatten()]
            f.write(",".join(row) + "\n")



if run_matrix:
    #run many simulations and calculate fc matrix and global strength for each, changing params for each one
    seeds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
    K_values = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75]
    Gamma_values = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75]
    n_regions = 80
    header = ["Seed", "K", "Gamma", "GlobalFC"]
    with open('fc_results.csv', 'w') as f:
        # Write header
        f.write(",".join(header) + "\n")

        # Loop through parameter combinations
        for seed in seeds:
            for K in K_values:
                for Gamma in Gamma_values:

                    #update random seeds and BOLD parameters
                    model.params.update({'seed': seed, 'bold_Gamma': Gamma, 'bold_K': K})
                    model.randomState = np.random.RandomState(seed)
                    if hasattr(model, 'noise_input'):
                       model.noise_input = []  # Clear noise buffers

                    model.run()
                    fc_matrix = np.corrcoef(model.BOLD.BOLD.T, rowvar=False)
                    global_fc = np.mean(np.abs(fc_matrix[np.triu_indices_from(fc_matrix, k=1)]))

                    # save row: sim#, K, Gamma, global_fc
                    row = [str(seed), f"{K:.4f}", f"{Gamma:.4f}", f"{global_fc:.4f}"]
                    f.write(",".join(row) + "\n")

                    print(f"Seed {seed}: K={K:.2f}, Γ={Gamma:.2f}, FC={global_fc:.4f}")



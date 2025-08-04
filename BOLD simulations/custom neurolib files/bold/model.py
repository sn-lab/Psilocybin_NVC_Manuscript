import numpy as np

from .timeIntegration import simulateBOLD


class BOLDModel:
    """
    Balloon-Windkessel BOLD simulator class.
    BOLD activity is downsampled to 0.5 Hz by default.

    BOLD simulation results are saved in t_BOLD, BOLD instance attributes.
    """

    def __init__(self, N, dt, normalize_input=False, normalize_max=50, rho=0.34, alpha=0.32, V0=0.02, k1=3.72, k2=0.53, k3=0.53, Gamma=0.41, K=0.65, Tau=0.98):
        self.N = N
        self.dt = dt  # dt of input activity in ms
        self.samplingRate_NDt = int(round(2000 / dt))  # downsample (0.5 Hz fMRI sampling rate)

        self.normalize_input = normalize_input
        self.normalize_max = normalize_max
        # return arrays
        self.t_BOLD = np.array([], dtype="f", ndmin=2)
        self.BOLD = np.array([], dtype="f", ndmin=2)
        self.all_Rates = np.array([], dtype="f", ndmin=2)
        self.BOLD_chunk = np.array([], dtype="f", ndmin=2)

        self.idxLastT = 0  # Index of the last computed t

        # initialize BOLD model variables
        self.X_BOLD = np.ones((N,))
        # Vasso dilatory signal
        self.F_BOLD = np.ones((N,))
        # Blood flow
        self.Q_BOLD = np.ones((N,))
        # Deoxyhemoglobin
        self.V_BOLD = np.ones((N,))
        # Blood volume

        #BOLD parameters (set by ALNModel)
        self.rho = rho
        self.alpha = alpha
        self.V0 = V0
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.Gamma = Gamma
        self.K = K
        self.Tau = Tau


    def run(self, activity, append=False):
        """Runs the Balloon-Windkessel BOLD simulation.

        Parameters:
            :param activity:     Neuronal firing rate in Hz
        
        :param activity: Neuronal firing rate in Hz
        :type activity: numpy.ndarray
        """

        # Compute the BOLD signal for the chunk
        BOLD_chunk, self.X_BOLD, self.F_BOLD, self.Q_BOLD, self.V_BOLD = simulateBOLD(
            activity,
            self.dt * 1e-3,
            10000 * np.ones((self.N,)),
            X=self.X_BOLD,
            F=self.F_BOLD,
            Q=self.Q_BOLD,
            V=self.V_BOLD,
            rho=self.rho,
            alpha=self.alpha,
            V0=self.V0,
            k1=self.k1,
            k2=self.k2,
            k3=self.k3,
            Gamma=self.Gamma,
            K=self.K,
            Tau=self.Tau
        )

        # downsample BOLD
        BOLD_resampled = BOLD_chunk[
            :, self.samplingRate_NDt - np.mod(self.idxLastT - 1, self.samplingRate_NDt) :: self.samplingRate_NDt
        ]
        t_new_idx = self.idxLastT + np.arange(activity.shape[1])
        t_BOLD_resampled = (
            t_new_idx[self.samplingRate_NDt - np.mod(self.idxLastT - 1, self.samplingRate_NDt) :: self.samplingRate_NDt]
            * self.dt
        )

        if self.BOLD.shape[1] == 0:
            # add new data
            self.t_BOLD = t_BOLD_resampled
            self.BOLD = BOLD_resampled
        elif append is True:
            # append new data to old data
            self.t_BOLD = np.hstack((self.t_BOLD, t_BOLD_resampled))
            self.BOLD = np.hstack((self.BOLD, BOLD_resampled))
        else:
            # overwrite old data
            self.t_BOLD = t_BOLD_resampled
            self.BOLD = BOLD_resampled

        self.BOLD_chunk = BOLD_resampled

        self.idxLastT = self.idxLastT + activity.shape[1]

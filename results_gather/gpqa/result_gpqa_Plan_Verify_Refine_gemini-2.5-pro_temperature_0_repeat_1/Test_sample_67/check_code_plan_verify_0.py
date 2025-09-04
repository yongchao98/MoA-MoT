import numpy as np
from scipy.linalg import expm
from numpy.fft import fft, fftfreq

def check_answer_correctness():
    """
    This function numerically simulates the spin precession problem to verify the
    oscillation frequency of the expectation value of the z-component of the magnetic moment.
    """
    # --- 1. Define constants and operators ---
    # We can use normalized units for the simulation where hbar=1.
    # The final frequency will be in terms of gamma and B.
    hbar = 1.0
    # Let's choose arbitrary values for gamma and B to perform the simulation.
    gamma = 1.5
    B = 2.0

    # The expected angular frequency from the answer C is omega = gamma * B
    omega_expected = gamma * B

    # Pauli matrices
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Spin operators (S = hbar/2 * sigma)
    Sy = (hbar / 2) * sigma_y
    Sz = (hbar / 2) * sigma_z

    # --- 2. Define the Hamiltonian and initial state ---
    # The magnetic field is along the +Y direction: B_vec = B * y_hat
    # The Hamiltonian is H = -gamma * B * Sy
    H = -gamma * B * Sy

    # The particle is initially in the spin-up state along the Z-axis.
    # In the standard Z-basis, |+z> = [1, 0]^T
    psi_0 = np.array([[1], [0]], dtype=complex)

    # --- 3. Simulate the time evolution ---
    # To find the frequency accurately using FFT, we need to simulate for several periods
    # and use a sufficient number of points.
    # The expected period is T = 2*pi / omega.
    period_expected = 2 * np.pi / omega_expected
    t_final = 10 * period_expected  # Simulate for 10 periods
    num_points = 2048               # Number of time steps (power of 2 is good for FFT)
    t_array = np.linspace(0, t_final, num_points)

    # Array to store the expectation value of mu_z
    exp_mu_z = np.zeros(num_points)

    for i, t in enumerate(t_array):
        # The time evolution operator is U(t) = exp(-iHt/hbar)
        U_t = expm(-1j * H * t / hbar)

        # The state at time t is |psi(t)> = U(t)|psi(0)>
        psi_t = U_t @ psi_0

        # The expectation value of Sz is <Sz>(t) = <psi(t)|Sz|psi(t)>
        exp_Sz_t = (psi_t.conj().T @ Sz @ psi_t).real

        # The expectation value of mu_z is gamma * <Sz>
        exp_mu_z[i] = gamma * exp_Sz_t[0, 0]

    # --- 4. Analyze the result to find the frequency using FFT ---
    # The signal is the time series of <mu_z>
    signal = exp_mu_z
    
    # The time step (sample spacing)
    time_step = t_array[1] - t_array[0]
    
    # Compute the FFT of the signal
    yf = fft(signal)
    # Compute the frequency bins. fftfreq gives frequency f, not angular frequency omega.
    xf = fftfreq(num_points, time_step)
    
    # We are interested in the dominant positive frequency.
    # Find the index of the peak in the positive frequency part of the spectrum.
    # We ignore the DC component at index 0.
    positive_freq_indices = np.where(xf > 0)
    peak_index = np.argmax(np.abs(yf[positive_freq_indices]))
    
    # The dominant frequency f (in Hz) from the simulation
    f_numerical = xf[positive_freq_indices][peak_index]
    
    # Convert to angular frequency omega (in rad/s)
    omega_numerical = 2 * np.pi * f_numerical

    # --- 5. Check correctness ---
    # Compare the numerically found frequency with the proposed answer's frequency.
    # The proposed answer is C, which corresponds to omega = gamma * B.
    if np.isclose(omega_numerical, omega_expected, rtol=1e-3):
        # The numerical simulation confirms the analytical result.
        # The analytical derivation in the provided answer is correct.
        return "Correct"
    else:
        return (f"The answer is incorrect. The proposed frequency is {omega_expected:.4f} (from gamma*B). "
                f"The numerically calculated frequency is {omega_numerical:.4f}.")

# Execute the check
result = check_answer_correctness()
# This will print "Correct" if the logic holds, or an error message otherwise.
# print(result)
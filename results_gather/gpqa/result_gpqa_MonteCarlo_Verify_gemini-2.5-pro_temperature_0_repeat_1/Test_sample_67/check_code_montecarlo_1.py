import numpy as np
from scipy.linalg import expm
from numpy.fft import fft, fftfreq

def check_spin_oscillation_frequency():
    """
    This function verifies the oscillation frequency of a spin-1/2 particle's
    z-component of magnetic moment when the magnetic field is suddenly switched
    from the Z-direction to the Y-direction.
    """
    # --- 1. Define Physical System and Constants ---
    # We can use normalized units (gamma=1, B=1, h_bar=1) for simplicity,
    # as the final frequency is proportional to gamma * B.
    gamma = 1.0
    B = 1.0
    h_bar = 1.0

    # Pauli matrices are the basis for spin-1/2 operators.
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Spin operators (S = h_bar/2 * sigma)
    Sy = (h_bar / 2) * sigma_y
    Sz = (h_bar / 2) * sigma_z

    # The new Hamiltonian is H = -mu . B = -gamma * S . B = -gamma * B * Sy
    H = -gamma * B * Sy

    # The initial state is spin-up along the Z-axis: |psi(0)> = |up_z>
    psi_0 = np.array([[1], [0]], dtype=complex)

    # --- 2. Analytical Derivation (Larmor Precession) ---
    # The spin precesses around the new B-field (Y-axis) at the Larmor frequency.
    # The expectation value <Sz(t)> oscillates as cos(omega*t).
    # The Larmor frequency is omega = gamma * B.
    expected_angular_frequency = gamma * B

    # --- 3. Numerical Simulation ---
    # We simulate the system's evolution to find the oscillation frequency numerically.
    
    # Set up time array for simulation. The period of one oscillation is T = 2*pi/omega.
    # We simulate for several periods to get a clean frequency peak.
    period = 2 * np.pi / expected_angular_frequency
    num_periods = 10
    num_points = 2048  # Use a power of 2 for FFT efficiency
    t_array = np.linspace(0, num_periods * period, num_points, endpoint=False)
    
    # Calculate the expectation value of Sz at each time step.
    sz_expectation_values = []
    for t in t_array:
        # The time evolution operator is U(t) = exp(-i*H*t/h_bar)
        U_t = expm(-1j * H * t / h_bar)
        
        # The state at time t is psi(t) = U(t) * psi(0)
        psi_t = U_t @ psi_0
        
        # The expectation value <Sz(t)> = <psi(t)|Sz|psi(t)>
        # np.vdot(a, b) computes a.conj().T @ b
        sz_exp_val = np.vdot(psi_t, Sz @ psi_t).real
        sz_expectation_values.append(sz_exp_val)

    # --- 4. Frequency Analysis using FFT ---
    # Perform a Fast Fourier Transform on the <Sz(t)> signal to find its frequency components.
    yf = fft(sz_expectation_values)
    
    # The fftfreq function gives the frequency bins corresponding to the FFT output.
    # The second argument is the time step (sample spacing).
    time_step = t_array[1] - t_array[0]
    xf = fftfreq(num_points, time_step)

    # We are interested in the dominant positive frequency.
    # Find the index of the peak in the positive frequency part of the spectrum.
    # We ignore the zero-frequency (DC) component at index 0.
    peak_index = np.argmax(np.abs(yf[1:num_points//2])) + 1
    
    # The frequency from the simulation is given in Hz (f).
    simulated_freq_hz = xf[peak_index]
    
    # Convert from frequency (f) to angular frequency (omega = 2*pi*f).
    simulated_angular_frequency = 2 * np.pi * simulated_freq_hz

    # --- 5. Verification ---
    # The provided answer is C, which corresponds to a frequency of gamma*B.
    answer_key = 'C'
    options = {
        'A': gamma * B / 2,
        'B': np.sqrt(2) * gamma * B,
        'C': gamma * B,
        'D': gamma * B / np.sqrt(2)
    }
    provided_answer_value = options[answer_key]

    # Check if the numerically simulated frequency matches the analytical result.
    if not np.isclose(simulated_angular_frequency, expected_angular_frequency, rtol=1e-3):
        return (f"The numerical simulation is inconsistent with the analytical derivation. "
                f"Analytical Frequency: {expected_angular_frequency:.4f}, "
                f"Simulated Frequency: {simulated_angular_frequency:.4f}")

    # Check if the provided answer matches the correct, verified frequency.
    if np.isclose(provided_answer_value, expected_angular_frequency):
        return "Correct"
    else:
        return (f"Incorrect. The provided answer '{answer_key}' corresponds to a frequency of "
                f"{provided_answer_value:.4f}, but the correct oscillation frequency, "
                f"derived analytically and verified by simulation, is {expected_angular_frequency:.4f} (gamma*B).")

# Run the check
result = check_spin_oscillation_frequency()
print(result)
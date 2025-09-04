import numpy as np
from scipy.linalg import expm
from scipy.fft import fft, fftfreq

def check_oscillation_frequency():
    """
    Numerically simulates the spin-1/2 particle's dynamics to verify the
    oscillation frequency of the z-component of its magnetic moment.
    """
    # --- 1. Define constants and parameters ---
    # We can set constants to simple values for the simulation, as the
    # relationship between them is what matters.
    hbar = 1.0
    # Let's choose arbitrary values for gamma and B
    gamma = 2.5
    B = 4.0
    
    # The answer D suggests the angular frequency is omega = gamma * B
    expected_omega = gamma * B
    
    # --- 2. Define spin operators (in the Z-basis) ---
    # S = (hbar/2) * sigma, where sigma are the Pauli matrices.
    Sx = (hbar / 2) * np.array([[0, 1], [1, 0]], dtype=complex)
    Sy = (hbar / 2) * np.array([[0, -1j], [1j, 0]], dtype=complex)
    Sz = (hbar / 2) * np.array([[1, 0], [0, -1]], dtype=complex)
    
    # --- 3. Define the initial state ---
    # The particle is initially aligned with a B-field in the +Z direction.
    # This means it's in the lowest energy eigenstate of H_initial = -gamma*B*Sz.
    # The eigenvalues are -gamma*B*(+/- hbar/2). The lowest energy corresponds to
    # the +hbar/2 eigenstate of Sz, which is the spin-up state |+z>.
    psi_initial = np.array([1, 0], dtype=complex) # |+z> = [1, 0]
    
    # --- 4. Define the new Hamiltonian (for t > 0) ---
    # The new magnetic field is B_f = B * j_hat (Y-direction).
    # The new Hamiltonian is H = -mu . B_f = -gamma * B * Sy.
    H_new = -gamma * B * Sy
    
    # --- 5. Simulate the time evolution and measure <mu_z> ---
    # We need to simulate for long enough to capture the oscillation.
    # The expected period is T = 2*pi / omega. Let's simulate for 5 periods.
    expected_period = 2 * np.pi / expected_omega
    t_max = 5 * expected_period
    num_steps = 2048 # Use a power of 2 for FFT efficiency
    t_values = np.linspace(0, t_max, num_steps, endpoint=False)
    
    # Array to store the results of <mu_z>(t)
    mu_z_values = np.zeros(num_steps)
    
    for i, t in enumerate(t_values):
        # Calculate the time evolution operator U(t) = exp(-i*H*t/hbar)
        U_t = expm(-1j * H_new * t / hbar)
        
        # Evolve the state: |psi(t)> = U(t)|psi(0)>
        psi_t = U_t @ psi_initial
        
        # Calculate the expectation value of Sz: <Sz> = <psi(t)|Sz|psi(t)>
        # The result must be real, so we take the real part to discard tiny
        # imaginary numerical errors.
        exp_val_Sz = np.real(psi_t.conj().T @ Sz @ psi_t)
        
        # Calculate the expectation value of mu_z
        mu_z_values[i] = gamma * exp_val_Sz
        
    # --- 6. Find the frequency using FFT ---
    # Perform the FFT on the signal
    yf = fft(mu_z_values)
    # Get the corresponding frequencies (f in Hz)
    # The sample spacing is dt = t_max / num_steps
    xf = fftfreq(num_steps, t_max / num_steps)
    
    # Find the frequency corresponding to the largest peak in the positive spectrum
    # We ignore the DC component at index 0
    positive_freq_indices = np.where(xf > 0)
    peak_index = np.argmax(np.abs(yf[positive_freq_indices]))
    
    # The measured frequency f (in Hz)
    measured_f = xf[positive_freq_indices][peak_index]
    # Convert to angular frequency omega (in rad/s)
    measured_omega = 2 * np.pi * measured_f
    
    # --- 7. Compare and return the result ---
    # Check if the measured angular frequency matches the expected one from answer D.
    # We use np.isclose for robust floating-point comparison.
    if np.isclose(measured_omega, expected_omega):
        return "Correct"
    else:
        return (f"The answer is incorrect. "
                f"The provided answer implies an angular frequency of omega = gamma * B = {expected_omega:.4f}. "
                f"However, the numerical simulation shows an angular frequency of omega = {measured_omega:.4f}. "
                f"The reasoning in the provided answer is correct, but the code check is designed to be an independent verification.")

# Run the check
result = check_oscillation_frequency()
print(result)
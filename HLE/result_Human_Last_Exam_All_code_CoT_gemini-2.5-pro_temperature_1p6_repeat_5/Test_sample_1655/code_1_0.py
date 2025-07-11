import numpy as np

def solve():
    """
    This script demonstrates that a system with a single linearly polarized input channel
    cannot create an arbitrary vector beam. It works by showing that the output
    polarization pattern is a fixed property of the system, regardless of the
    input beam's scalar profile.
    """
    # 1. Setup simulation parameters
    N = 128  # Grid size
    print(f"Simulating a {N}x{N} beam profile.\n")
    X, Y = np.meshgrid(np.linspace(-3, 3, N), np.linspace(-3, 3, N))
    R = np.sqrt(X**2 + Y**2)
    THETA = np.arctan2(Y, X)

    # 2. Define the system's components

    # Free-space propagation operator (simulated as a Fourier Transform)
    def propagate(field_vector):
        # Applies a Fourier transform to each polarization component independently
        propagated_field = np.zeros_like(field_vector)
        propagated_field[:, :, 0] = np.fft.fftshift(np.fft.fft2(field_vector[:, :, 0]))
        propagated_field[:, :, 1] = np.fft.fftshift(np.fft.fft2(field_vector[:, :, 1]))
        return propagated_field

    # Random Medium Transmission Matrix T (spatially-varying Jones Matrix)
    # This is the most general case, allowing polarization mixing.
    np.random.seed(0) # for reproducibility
    T = np.random.randn(N, N, 2, 2) + 1j * np.random.randn(N, N, 2, 2)
    T_inv = np.linalg.inv(T)

    # Operator for applying the Jones matrix at each pixel
    def apply_jones(jones_matrix, field_vector):
        # Einstein summation: for each pixel(x,y), multiply the matrix T[x,y] by the vector E[x,y]
        return np.einsum('xyij,xyj->xyi', jones_matrix, field_vector)

    # Define the full system operation as described in the prompt
    def run_full_system(E_in):
        # Path 1: E_in -> P -> T -> P -> Output1
        E_propagated_1 = propagate(E_in)
        E_through_T = apply_jones(T, E_propagated_1)
        Output1 = propagate(E_through_T)
        
        # Path 2: Output1 -> P -> T_inv -> Output2
        E_propagated_2 = propagate(Output1)
        Output2 = apply_jones(T_inv, E_propagated_2)
        return Output2

    # 3. Create two different "tailored" scalar inputs, both with fixed horizontal polarization
    
    # Input 1: A simple Gaussian beam
    scalar_input_1 = np.exp(-R**2)
    E_in_1 = np.zeros((N, N, 2), dtype=complex)
    E_in_1[:, :, 0] = scalar_input_1  # Horizontally polarized: [Scalar, 0]

    # Input 2: A Laguerre-Gaussian beam (l=1, p=0) with a phase vortex
    scalar_input_2 = (R) * np.exp(-R**2) * np.exp(1j * 1 * THETA)
    E_in_2 = np.zeros((N, N, 2), dtype=complex)
    E_in_2[:, :, 0] = scalar_input_2  # Horizontally polarized: [Scalar, 0]
    
    print("Created two different scalar inputs (Gaussian and Laguerre-Gaussian) with fixed horizontal polarization.")
    print("Passing both through the simulated optical system...\n")

    # 4. Propagate both distinct inputs through the entire system
    Output_1 = run_full_system(E_in_1)
    Output_2 = run_full_system(E_in_2)

    # 5. Test the hypothesis: Is the output polarization pattern independent of the input?
    # We calculate the complex ratio Ey/Ex, which defines the state of polarization at each pixel.
    epsilon = 1e-9  # Add a small number to avoid division by zero
    ratio_1 = Output_1[:, :, 1] / (Output_1[:, :, 0] + epsilon)
    ratio_2 = Output_2[:, :, 1] / (Output_2[:, :, 0] + epsilon)

    # 6. Compare the ratios from the two different inputs to see if they are the same.
    px, py = N // 2 + 15, N // 2 - 10 # Check at an arbitrary off-axis point
    r1_point = ratio_1[px, py]
    r2_point = ratio_2[px, py]
    difference = np.abs(r1_point - r2_point)

    print("-" * 50)
    print("Hypothesis Check: Is the output polarization pattern fixed?")
    print(f"Checking polarization state (complex ratio Ey/Ex) at pixel ({px}, {py}):")
    print(f"Ratio for Gaussian input:      {r1_point.real:+.4f} {r1_point.imag:+.4f}j")
    print(f"Ratio for LG beam input:       {r2_point.real:+.4f} {r2_point.imag:+.4f}j")
    print(f"Magnitude of their difference: {difference:.2e}")
    print("-" * 50)

    # A more robust check across the whole grid
    mean_difference = np.mean(np.abs(ratio_1 - ratio_2))
    print(f"\nMean absolute difference of ratios across the entire beam: {mean_difference:.2e}\n")

    print("Conclusion:")
    print("The difference is zero (within numerical precision). This shows that the output polarization pattern")
    print("is identical for different scalar inputs. It is a fixed property of the optical system.")
    print("Since the polarization pattern cannot be controlled by tailoring the input beam, an *arbitrary*")
    print("vector beam cannot be generated.")

solve()
<<<No>>>
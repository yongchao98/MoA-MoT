import cmath

def calculate_bare_greens_function():
    """
    Calculates the bare Green's function G_0 for given parameters
    and prints the functional dependence and final equation.
    """
    # Define parameters for the calculation
    # (using natural units where h-bar = 1)
    omega = 3.5  # Frequency (energy of the probe)
    epsilon_k = 3.0  # Single-particle energy eigenvalue
    eta = 0.1  # Infinitesimal positive value for causality

    # The formula for the bare Green's function
    # G_0(k, omega) = 1 / (omega - epsilon_k + i*eta)
    
    # In Python, the imaginary unit 'i' is represented by 'j'
    denominator = (omega - epsilon_k) + 1j * eta
    g0 = 1 / denominator

    # Output the explanation and the final equation with numbers
    print("The functional dependence of the bare Green's function G_0 on the single-particle energy epsilon_k is:")
    print("G_0(k, omega) = 1 / (omega - epsilon_k + i*eta)\n")
    
    print("This shows G_0 is inversely proportional to the difference between the frequency omega and the energy epsilon_k.")
    print("\n--- Example Calculation ---")
    print(f"Let's use the following values:")
    print(f"omega (frequency) = {omega}")
    print(f"epsilon_k (single-particle energy) = {epsilon_k}")
    print(f"eta (infinitesimal for causality) = {eta}\n")

    # Print the final equation with the numbers substituted in
    print("The equation is:")
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + i * {eta})")
    
    # Simplify the denominator to show the next step
    real_part_denom = omega - epsilon_k
    print(f"G_0 = 1 / ({real_part_denom:.2f} + {eta:.2f}i)\n")
    
    # Print the final complex number result
    print(f"The resulting value for the bare Green's function is:")
    print(f"G_0 = {g0.real:.4f} + ({g0.imag:.4f})i")

# Run the calculation
calculate_bare_greens_function()
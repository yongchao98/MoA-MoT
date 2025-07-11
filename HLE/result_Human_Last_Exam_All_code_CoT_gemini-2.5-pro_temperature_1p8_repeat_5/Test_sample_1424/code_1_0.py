import cmath

def calculate_bare_green_function():
    """
    Explains and calculates the bare Green's function G_0.

    The bare Green's function G_0(k, omega) describes the propagation of a
    non-interacting particle. Its functional dependence on the single-particle
    energy eigenvalue epsilon_k is given by:

    G_0(k, omega) = 1 / (omega - epsilon_k + i * eta * sign(epsilon_k - mu))
    """

    # --- Parameters for a numerical example ---
    # We choose some arbitrary but physically plausible values for demonstration.
    omega = 2.5      # Probe energy (e.g., in eV)
    epsilon_k = 2.0  # Single-particle energy eigenvalue (e.g., in eV)
    mu = 0.0         # Chemical potential (we assume epsilon_k > mu, i.e., an unoccupied state)
    eta = 0.1        # A small positive infinitesimal for causality

    # --- Explanation ---
    print("The functional dependence of the bare Green's function G_0 on the single-particle energy epsilon_k is:")
    print("G_0(k, omega) = 1 / (omega - epsilon_k + i * eta * sgn(epsilon_k - mu))\n")
    print("Let's calculate this for a specific example:")
    print(f"  Probe Energy (omega)      = {omega}")
    print(f"  Single-particle (epsilon_k) = {epsilon_k}")
    print(f"  Chemical Potential (mu)   = {mu}")
    print(f"  Infinitesimal (eta)       = {eta}\n")

    # --- Step-by-step calculation ---

    # 1. Determine the sign term
    if (epsilon_k - mu) > 0:
        sign_val = 1
    elif (epsilon_k - mu) < 0:
        sign_val = -1
    else:
        # This case is technically singular without eta, but we can assign 0 for the calculation
        sign_val = 0
    print(f"Step 1: Calculate the sign term sgn(epsilon_k - mu)")
    print(f"   sgn({epsilon_k} - {mu}) = sgn({epsilon_k - mu}) = {sign_val}\n")

    # 2. Calculate the terms in the denominator
    real_part = omega - epsilon_k
    imag_part = eta * sign_val
    print("Step 2: Calculate the real and imaginary parts of the denominator")
    print(f"   Real Part: omega - epsilon_k = {omega} - {epsilon_k} = {real_part}")
    print(f"   Imaginary Part: eta * sgn(...) = {eta} * {sign_val} = {imag_part}\n")

    # 3. Construct the denominator and the full equation
    denominator = complex(real_part, imag_part)
    print("Step 3: Assemble the final equation")
    print(f"   G_0 = 1 / (({real_part}) + ({imag_part})i)")
    print(f"   G_0 = 1 / {denominator}\n")

    # 4. Calculate the final result by taking the inverse
    G_0 = 1 / denominator
    print("Step 4: Calculate the final value of the Green's function")
    print(f"   G_0 = {G_0.real:.4f} + ({G_0.imag:.4f})i")


# Run the calculation and print the results
calculate_bare_green_function()
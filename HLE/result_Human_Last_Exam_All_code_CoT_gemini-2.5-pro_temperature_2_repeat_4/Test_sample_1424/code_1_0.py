import cmath

def calculate_bare_green_function(omega, epsilon_k, eta):
    """
    Calculates the bare Green's function G_0(k, ω).

    This function demonstrates the dependence of G_0 on the single-particle
    energy eigenvalue ϵ_k.

    Args:
        omega (float): The frequency (energy) of the particle.
        epsilon_k (float): The single-particle energy eigenvalue.
        eta (float): A small positive infinitesimal for broadening/causality.

    Returns:
        tuple: A tuple containing the complex value of G_0 and the denominator.
    """
    # The denominator of the Green's function formula, which shows the
    # direct dependence on the difference between ω and ϵ_k.
    denominator = (omega - epsilon_k) + 1j * eta
    
    # The Green's function is the inverse of the denominator.
    G0 = 1 / denominator
    return G0, denominator

# --- Example Parameters ---
# Let's assume energy units are in eV (electron-volts)
omega = 2.0      # An example frequency (energy of observation)
epsilon_k = 1.5  # The energy of the single-particle state k
eta = 0.1        # A small positive constant representing lifetime effects or for calculation purposes

# --- Calculation ---
G0_result, denominator_result = calculate_bare_green_function(omega, epsilon_k, eta)

# --- Output the Explanation and Final Equation ---
print("This script calculates the bare Green's function G_0.")
print("The functional dependence of G_0 on the single-particle energy ϵ_k is given by the formula:")
print("G_0(ω, ϵ_k) = 1 / (ω - ϵ_k + iη)")
print("-" * 60)

# Print each number used in the final equation
print("Given the parameters:")
print(f"  Frequency ω          = {omega}")
print(f"  Single-particle ϵ_k  = {epsilon_k}")
print(f"  Broadening factor η  = {eta}\n")

print("Step 1: Calculate the denominator of the equation.")
# We use .real and .imag to explicitly show the numbers that form the complex denominator
print(f"Denominator = (ω - ϵ_k) + iη")
print(f"Denominator = ({omega} - {epsilon_k}) + {eta}j")
print(f"Denominator = {denominator_result.real} + {denominator_result.imag}j\n")

print("Step 2: Calculate G_0 by taking the inverse of the denominator.")
print(f"G_0 = 1 / ({denominator_result})")
print(f"Final Result: G_0 = {G0_result}")
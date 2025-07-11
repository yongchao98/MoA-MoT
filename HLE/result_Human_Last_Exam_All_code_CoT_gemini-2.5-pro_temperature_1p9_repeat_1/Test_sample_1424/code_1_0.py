import cmath

def calculate_bare_greens_function(energy_probe, energy_eigenvalue, delta):
    """
    Calculates the bare Green's function G_0 in the energy domain.

    The formula is G_0 = 1 / (E - epsilon_k + i*delta)

    Args:
        energy_probe (float): The energy 'E' at which the system is probed.
        energy_eigenvalue (float): The single-particle energy eigenvalue 'epsilon_k'.
        delta (float): A small positive infinitesimal.

    Returns:
        complex: The value of the bare Green's function G_0.
    """
    # In Python, the imaginary unit 'i' is represented as '1j'
    denominator = (energy_probe - energy_eigenvalue) + 1j * delta
    g0 = 1 / denominator
    return g0

# --- Parameters for the calculation ---
# Let's assume the particle has a single-particle energy eigenvalue of 2.0 eV.
epsilon_k = 2.0

# Let's probe the system at an energy of 2.1 eV, slightly off resonance.
E = 2.1

# We use a small positive value for delta.
delta = 0.05

# --- Perform the calculation ---
g0_result = calculate_bare_greens_function(E, epsilon_k, delta)

# --- Print the results ---
print("The functional dependence of the bare Green's function G_0 on the single-particle energy eigenvalue epsilon_k is:")
print("G_0(E) = 1 / (E - epsilon_k + i*delta)\n")

print("Calculating for the given values:")
print(f"E = {E}")
print(f"epsilon_k = {epsilon_k}")
print(f"delta = {delta}\n")

# Display the final equation with numbers substituted
# The '.2f' formats the numbers to two decimal places for clarity.
print("Final Equation:")
print(f"G_0 = 1 / ({E:.2f} - {epsilon_k:.2f} + i*{delta:.2f})")
print(f"G_0 = 1 / ({(E - epsilon_k):.2f} + {delta:.2f}i)")


# Print the final numerical answer
# We print the real and imaginary parts for clarity.
print(f"\nResult: G_0 = {g0_result:.4f}")

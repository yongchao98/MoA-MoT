import cmath

def calculate_bare_greens_function(E, epsilon_k, eta=0.1):
    """
    Calculates the retarded bare Green's function G_0(E, epsilon_k).

    Args:
        E (float): The energy of the probing particle.
        epsilon_k (float): The single-particle energy eigenvalue for state k.
        eta (float): A small positive constant for causality.

    Returns:
        complex: The value of the bare Green's function.
    """
    denominator = (E - epsilon_k) + 1j * eta
    G0 = 1 / denominator
    return G0

# --- Main execution ---

# Define the single-particle energy eigenvalue (in arbitrary units, e.g., eV)
epsilon_k = 2.0

# Define the probe energy, set close to the particle energy to see a strong response
E = 2.05

# Define the infinitesimal constant eta (kept small but finite for calculation)
eta = 0.1

# Calculate the Green's function
G0 = calculate_bare_greens_function(E, epsilon_k, eta)

# Print the explanation and the result
print("The retarded bare Green's function G_0 depends on the probe energy (E) and the single-particle energy (epsilon_k).")
print("The formula is: G_0(k, E) = 1 / (E - epsilon_k + i*eta)")
print("\n--- Calculation Example ---")
print(f"Given values:")
print(f"  Probe Energy (E)          = {E}")
print(f"  Particle Energy (epsilon_k) = {epsilon_k}")
print(f"  Constant (eta)            = {eta}")
print("\nSubstituting these values into the formula:")
# Here we output each number in the final equation
print(f"G_0 = 1 / ({E} - {epsilon_k} + i*{eta})")
print(f"G_0 = 1 / ({(E - epsilon_k):.2f} + {eta:.2f}i)")
print(f"\nThe calculated complex value is:")
print(f"G_0 = {G0.real:.4f} + {G0.imag:.4f}i")

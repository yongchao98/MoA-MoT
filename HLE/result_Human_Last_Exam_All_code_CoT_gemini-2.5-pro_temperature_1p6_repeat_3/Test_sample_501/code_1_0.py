import math

def calculate_adiabatic_polymer_force(x, l, n, E0):
    """
    Calculates the attractive force between the ends of a thermally isolated
    freely jointed polymer chain.

    Args:
        x (float): The end-to-end separation of the polymer.
        l (float): The length of a single strut (link).
        n (int): The number of links in the chain.
        E0 (float): The kinetic energy of the polymer at zero extension (x=0).

    Returns:
        float: The attractive force F(x). A negative value indicates a
               restoring force pulling the ends together.
    """
    if n <= 0:
        raise ValueError("Number of links n must be positive.")
    if l <= 0:
        raise ValueError("Link length l must be positive.")

    # The force law is F(x) = - (3 * E0 * x / (n^2 * l^2)) * exp((3 * x^2) / (2 * n^2 * l^2))
    
    # Calculate the pre-factor (linear term)
    prefactor = (3 * E0 * x) / (n**2 * l**2)
    
    # Calculate the term in the exponent
    exponent_term = (3 * x**2) / (2 * n**2 * l**2)
    
    # Calculate the full force
    force = -prefactor * math.exp(exponent_term)
    
    return force, prefactor, exponent_term

# --- Example Calculation ---

# Define the polymer parameters
n_links = 1000  # Number of links
link_length = 1.0e-9  # meters (1 nm)

# Define the initial energy and extension
# We can estimate E(0) from an initial temperature T0, e.g., 300 K
# E(0) = n * k_B * T0, where k_B is the Boltzmann constant
boltzmann_constant = 1.380649e-23  # J/K
initial_temp = 300  # Kelvin
E_initial = n_links * boltzmann_constant * initial_temp # Joules

# Set the extension. This should be small compared to the total contour length (n*l)
# Total length L = 1000 * 1e-9 = 1e-6 m.
extension_x = 10e-9 # meters (10 nm)

# Calculate the force
force, prefactor, exponent_term = calculate_adiabatic_polymer_force(extension_x, link_length, n_links, E_initial)

# --- Output Results ---
print("--- Polymer Parameters ---")
print(f"Number of links (n): {n_links}")
print(f"Length of each link (l): {link_length:.2e} m")
print(f"Initial kinetic energy (E(0)): {E_initial:.3e} J (corresponds to T={initial_temp} K)")
print(f"End-to-end separation (x): {extension_x:.2e} m")
print("\n--- Force Calculation ---")
print("Force Law: F(x) = - (3*E(0)*x / (n^2*l^2)) * exp(3*x^2 / (2*n^2*l^2))")
print("\n--- Calculated Values ---")
print(f"Linear term value (3*E(0)*x / (n^2*l^2)): {prefactor:.3e} N")
print(f"Exponent term value (3*x^2 / (2*n^2*l^2)): {exponent_term:.3e}")
print(f"Final attractive force F(x): {force:.3e} N")

# The final answer is the derived force law formula
final_answer = "F(x) = - (3 * E(0) * x / (n**2 * l**2)) * exp((3 * x**2) / (2 * n**2 * l**2))"
# The provided code implements this formula. No single numerical answer is possible
# as it's a functional relationship. For the example values given, the force is the calculated value.
# We will present the symbolic formula as the final answer.
print(f"\n<<<F(x) = -(3*E(0)*x/(n^2*l^2)) * exp(3*x^2/(2*n^2*l^2))>>>")
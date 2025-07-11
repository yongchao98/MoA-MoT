import math
import scipy.constants

def calculate_joule_heat():
    """
    Calculates the Joule heat dissipated by a leaking and shrinking charged sphere.

    The function assumes the charge leakage is the dominant process, meaning the
    total dissipated heat is equal to the initial electrostatic potential energy
    of the sphere.
    """
    # --- User-defined variables ---
    # You can change these values to match your specific problem.
    # Radius of the sphere in meters.
    a = 0.1
    # Initial potential of the sphere in volts.
    V = 10000.0

    # --- Physical constant ---
    # Permittivity of free space in Farads per meter (F/m).
    epsilon_0 = scipy.constants.epsilon_0

    # --- Calculations ---
    # The formula for the total heat dissipated (H) is derived from the initial
    # electrostatic energy (U) of the sphere: H = U = 2 * pi * epsilon_0 * a * V^2.

    # 1. Initial capacitance of the sphere (C = 4 * pi * epsilon_0 * a)
    capacitance = 4 * math.pi * epsilon_0 * a

    # 2. Initial stored electrostatic energy (U = 1/2 * C * V^2)
    # This energy is fully converted to Joule heat in our model.
    heat_dissipated = 0.5 * capacitance * V**2
    
    # An alternative direct calculation using the final formula
    # heat_dissipated_direct = 2 * math.pi * epsilon_0 * a * V**2

    # --- Output ---
    print("This script calculates the Joule heat dissipated from a charged sphere.")
    print("It assumes the heat (H) equals the initial electrostatic energy (U_initial).\n")
    print("Given values:")
    print(f"  Initial radius (a) = {a} m")
    print(f"  Initial potential (V) = {V} V")
    print(f"  Permittivity of free space (ε₀) = {epsilon_0:.4e} F/m\n")

    print("Step-by-step calculation:")
    # Using the equation H = 2 * pi * ε₀ * a * V²
    # We output each number in the final equation as requested.
    print(f"H = 2 * π * ε₀ * a * V²")
    print(f"H = 2 * {math.pi:.5f} * {epsilon_0:.5e} * {a} * {V**2:.2e}")
    print(f"H = {heat_dissipated:.5e} Joules")

# Execute the function
calculate_joule_heat()
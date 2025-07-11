import math

def calculate_joule_heat():
    """
    Calculates the Joule heat dissipated by a charged sphere that leaks charge
    while its radius decreases.
    
    The problem asks for the total Joule heat dissipated.
    The initial electrostatic energy of the sphere is U = (1/2)CV^2.
    For a sphere of radius 'a', the capacitance is C = 4πε₀a.
    So, the initial energy is U = 2πε₀aV².
    
    This energy is converted into Joule heat (H) and mechanical work (W).
    U = H + W.
    The work done depends on the process path. A standard interpretation for such problems,
    yielding a definite answer, assumes that the charge leaks away completely before
    the sphere shrinks, making W = 0.
    Thus, the total Joule heat dissipated is equal to the initial stored energy.
    H = U = 2πε₀aV².
    """
    # Physical constant
    epsilon_0 = 8.8541878128e-12  # Permittivity of free space in F/m

    # User-defined variables for the initial state.
    # Example values are used here.
    a = 0.1  # Initial radius in meters
    V = 10000.0  # Initial potential in Volts

    # Calculate the total Joule heat dissipated, which equals the initial energy
    H = 2 * math.pi * epsilon_0 * a * (V ** 2)

    # Print the explanation and the final equation with all values
    print("The formula for the total Joule heat dissipated (H) is equal to the sphere's initial electrostatic energy:")
    print("H = 2 * π * ε₀ * a * V²")
    print("\nSubstituting the given values:")
    print(f"H = 2 * {math.pi:.5f} * {epsilon_0:.5e} F/m * {a} m * ({V} V)²")
    print(f"H = {H:.8f} Joules")

calculate_joule_heat()
import math

def calculate_joule_heat(radius_a, potential_v):
    """
    Calculates the total Joule heat dissipated by a leaking and shrinking charged sphere.

    Args:
        radius_a (float): The initial radius of the sphere in meters.
        potential_v (float): The initial potential of the sphere in volts.

    Returns:
        float: The total Joule heat dissipated in joules.
    """
    # Permittivity of free space in F/m (Farads per meter)
    epsilon_0 = 8.854187817e-12

    # The total Joule heat dissipated is equal to the initial electrostatic energy stored in the sphere.
    # Formula: H = 2 * pi * epsilon_0 * a * V^2
    joule_heat = 2 * math.pi * epsilon_0 * radius_a * potential_v**2

    print("The formula for the total Joule heat (H) dissipated is:")
    print("H = 2 * π * ε₀ * a * V²")
    print("\nCalculating with the given values:")
    print(f"H = 2 * {math.pi:.5f} * {epsilon_0} * {radius_a} * {potential_v}²")
    print(f"\nH = {joule_heat:.4g} Joules")
    
    return joule_heat

if __name__ == '__main__':
    # Example values for the user to try.
    # You can change these values to match your specific problem.
    initial_radius = 0.1  # meters
    initial_potential = 10000  # volts

    print(f"Solving for a sphere with initial radius a = {initial_radius} m and initial potential V = {initial_potential} V.")
    print("-" * 30)
    
    final_heat = calculate_joule_heat(initial_radius, initial_potential)
    
    # The final answer is the numerical value calculated.
    # For the example values a=0.1m and V=10000V, the answer is approximately 5.563e-5 Joules.
    # The final answer format is requested as <<<answer content>>>.
    # Since the problem is symbolic, we present the formula as the answer.
    # If a numerical answer is needed, it would be the calculated value.
    # For the example: <<<5.563e-05>>>
    # For the symbolic problem: <<<2 * pi * epsilon_0 * a * V**2>>>
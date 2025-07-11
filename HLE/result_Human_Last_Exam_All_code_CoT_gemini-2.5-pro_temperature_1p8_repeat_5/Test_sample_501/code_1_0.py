import math

def calculate_polymer_force(x, n, l, E0):
    """
    Calculates the restorative force of a thermally isolated freely jointed chain.

    Args:
        x (float): The end-to-end separation of the polymer.
        n (int): The number of mass points (or segments) in the polymer chain.
        l (float): The length of a single strut between mass points.
        E0 (float): The kinetic energy of the polymer at zero extension (x=0).
    
    Returns:
        float: The force F(x).
    """
    if n <= 0 or l <= 0:
        raise ValueError("n and l must be positive.")

    # Calculate the pre-factor
    prefactor = (2 * E0 * x) / (n**2 * l**2)
    
    # Calculate the exponential term
    exponent = - (x**2) / (n**2 * l**2)
    exp_term = math.exp(exponent)
    
    # The force is restorative, hence the negative sign
    force = -prefactor * exp_term
    
    return force

def main():
    """
    Main function to demonstrate the calculation with example values.
    """
    # Example parameters
    n = 1000  # Number of segments (large)
    l = 1.0   # Length of one segment (e.g., in nm)
    x = 50.0  # End-to-end separation (e.g., in nm), small compared to n*l
    E0 = 4.14e-21 * 100 # Initial kinetic energy (e.g., 100 * k_B * 300K, in Joules)

    print(f"Calculating the force for a polymer with the following parameters:")
    print(f"Number of segments (n): {n}")
    print(f"Length of a segment (l): {l}")
    print(f"End-to-end separation (x): {x}")
    print(f"Initial kinetic energy (E(0)): {E0:.4g}")
    print("-" * 30)

    # Calculate the force
    force = calculate_polymer_force(x, n, l, E0)

    # Print the equation with values plugged in
    print("Force Law: F(x) = - (2 * E(0) * x / (n^2 * l^2)) * exp(-x^2 / (n^2 * l^2))")
    print("\nPlugging in the values:")
    
    # Breaking down the calculation for clarity
    prefactor_val = (2 * E0 * x) / (n**2 * l**2)
    exponent_val = - (x**2) / (n**2 * l**2)

    print(f"F({x}) = - (2 * {E0:.4g} * {x}) / ({n}^2 * {l}^2) * exp(-{x}^2 / ({n}^2 * {l}^2))")
    print(f"F({x}) = - ({prefactor_val:.4g}) * exp({exponent_val:.4g})")
    print(f"F({x}) = {force:.4g} Newtons")

if __name__ == "__main__":
    main()
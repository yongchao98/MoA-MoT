import math

def evaluate_nuclear_attraction(zeta_i, zeta_j):
    """
    Evaluates the integral <phi_i | 1/r | phi_j> for 1s Slater orbitals.

    Args:
        zeta_i (float): Orbital exponent for phi_i.
        zeta_j (float): Orbital exponent for phi_j.

    Returns:
        float: The value of the integral.
    """
    # The integral I = 4 * (zeta_i * zeta_j)^(3/2) / (zeta_i + zeta_j)^2
    numerator = 4 * (zeta_i * zeta_j)**1.5
    denominator = (zeta_i + zeta_j)**2
    
    if denominator == 0:
        return float('inf') # Avoid division by zero, though unlikely with physical zeta values
        
    return numerator / denominator

def main():
    """
    Main function to set parameters, calculate the integral, and print results.
    """
    # Define the orbital exponents for two 1s Slater orbitals
    zeta_i = 1.0
    zeta_j = 1.2

    # Calculate the integral
    result = evaluate_nuclear_attraction(zeta_i, zeta_j)

    # Print the result showing the numbers used in the calculation
    # as per the derived formula: I = 4 * (z_i * z_j)^(3/2) / (z_i + z_j)^2
    print(f"For orbital exponents zeta_i = {zeta_i} and zeta_j = {zeta_j}:")
    print("The integral <phi_i | 1/r | phi_j> is evaluated as:")
    
    # We show the final equation with the numbers plugged in
    equation = f"I = (4 * ({zeta_i} * {zeta_j})**(3/2)) / ({zeta_i} + {zeta_j})**2"
    
    print(equation)
    print(f"I = {result}")

if __name__ == "__main__":
    main()

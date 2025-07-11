import math

def print_solution():
    """
    This function prints the derived analytical solution for the shock wave density profile.
    """

    # Define the numerical constants present in the final equation.
    c1 = 1
    c2 = 4
    c3 = 3
    c4 = 8

    # Provide context for the solution
    print("The analytical solution for the density profile of a Mach sqrt(3) shock wave in a monatomic ideal gas with Pr = 3/4 is derived from the Navier-Stokes equations.")
    print("The density (rho) is expressed in units of the ambient density (rho_0).")
    print("The position (x) is expressed in units of the ambient conductive length scale L, where L = kappa/(rho_0 * M * c_0 * C_v).")
    print("\nThe final equation for the normalized density profile is:")
    print("-" * 50)

    # Construct and print the final equation as a formatted string.
    equation_str = f"rho(x)/rho_0 = {c1} + sqrt(exp({c2}/{c3} * x/L) / (exp({c2}/{c3} * x/L) + {c4}))"
    print(equation_str)
    
    print("-" * 50)

    # As requested, output each number in the final equation.
    print("\nThe numerical constants in the equation are:")
    print(f"1. The additive term outside the square root: {c1}")
    print(f"2. The coefficient in the exponential term is a fraction: {c2}/{c3}")
    print(f"3. The constant added to the exponential in the denominator of the square root term: {c4}")

if __name__ == '__main__':
    print_solution()
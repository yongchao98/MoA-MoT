import numpy as np

def calculate_density_profile(xi):
    """
    Calculates the non-dimensional density profile of a Mach sqrt(3) shock wave.

    The solution is for a monatomic ideal gas (gamma = 5/3) with a constant
    dynamic viscosity and a Prandtl number of 3/4.

    Args:
        xi (float or np.ndarray): The non-dimensional position, xi = x / L,
            where L is the ambient conductive length scale.

    Returns:
        float or np.ndarray: The non-dimensional density rho' = rho / rho_0.
    """
    
    # For a Mach sqrt(3) shock with gamma=5/3 and Pr=3/4, the analytical
    # solution for the density profile rho'(xi) can be found explicitly.
    # The derivation leads to the following quadratic equation for rho':
    # rho'^2 * (3 + exp(k*xi)) - rho' * 2 * (3 + exp(k*xi)) + 3 = 0
    # where k = 5/3.
    # Solving this quadratic equation for rho' gives the explicit solution below.

    # The exponent in the solution is k*xi where k=5/3
    k = 5/3
    exponent_term = np.exp(k * xi)
    
    # The final explicit solution is rho'(xi) = 1 + sqrt( E / (3 + E) )
    # where E is the exponential term. We take the positive root to ensure
    # the density increases from 1 to 2 across the shock.
    density_profile = 1 + np.sqrt(exponent_term / (3 + exponent_term))
    
    return density_profile

# The problem asks for the analytical solution. The derived equation is:
# rho'(xi) = 1 + sqrt(exp((5/3) * xi) / (3 + exp((5/3) * xi)))
# The following code demonstrates how to use this solution.

print("The analytical solution for the non-dimensional density profile rho'(xi) is:")
print("rho'(xi) = 1 + (exp((5/3) * xi) / (3 + exp((5/3) * xi)))**0.5")
print("-" * 60)
print("Here, rho' = rho/rho_0 is the density normalized by the ambient density.")
print("And xi = x/L is the position normalized by the ambient conductive length scale.")
print("\nThe numbers in this equation are 1, 5, 3, 3, and 0.5 (from the square root).")
print("\nBelow are some example values calculated using a function that implements this solution:")

# Calculate and print the density for a few representative positions
# The origin xi=0 is defined as the point where rho' = 1.5
positions = np.array([-5.0, -1.0, 0.0, 1.0, 5.0])
print("\n  xi\t   |   rho'(xi)")
print("------------------------")
for xi_val in positions:
    rho_val = calculate_density_profile(xi_val)
    print(f" {xi_val: 6.2f}\t   |\t {rho_val:.6f}")
    
# We can also check the upstream and downstream boundary conditions
print("\nChecking the boundary conditions:")
print(f"As xi -> -infinity, rho' approaches {calculate_density_profile(-100):.6f} (The theoretical value is 1).")
print(f"As xi -> +infinity, rho' approaches {calculate_density_profile(100):.6f} (The theoretical value is 2).")

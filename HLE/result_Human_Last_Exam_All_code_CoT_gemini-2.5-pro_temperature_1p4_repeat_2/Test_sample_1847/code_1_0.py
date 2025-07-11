import math

def solve_density_profile():
    """
    Calculates the coefficients for the analytical solution of a Mach sqrt(3) shock wave
    in a monatomic ideal gas with Pr = 3/4.
    """

    # 1. Define physical parameters
    # For a monatomic ideal gas, the ratio of specific heats (gamma) is 5/3.
    gamma = 5.0 / 3.0
    # The shock Mach number squared is given as 3.
    M_sq = 3.0
    
    # 2. Calculate the post-shock to pre-shock density ratio (rho_1/rho_0)
    # This comes from the Rankine-Hugoniot jump conditions for a normal shock.
    rho_ratio = (gamma + 1) * M_sq / (2.0 + (gamma - 1) * M_sq)

    # 3. Determine the coefficients for the tanh solution
    # The solution for the normalized density profile rho(x)/rho_0 has the form:
    # rho(x)/rho_0 = A - B * tanh(C * x/L)
    # where L is the given characteristic length scale.
    
    # The coefficients A and B are determined by the boundary conditions:
    # At x -> -infinity, rho/rho_0 -> 1. Since tanh(-inf) = -1, we have A + B = 1.
    # At x -> +infinity, rho/rho_0 -> rho_ratio. Since tanh(+inf) = 1, we have A - B = rho_ratio.
    # Solving this system of two linear equations:
    A = (rho_ratio + 1.0) / 2.0
    B = (1.0 - rho_ratio) / 2.0  # This will be negative, matching A - B*tanh
    
    # 4. Determine the coefficient C in the argument of tanh
    # Based on the approximated shock-structure theory for Pr = 3/4, the coefficient is
    # C = 3 * (gamma + 1) / 16.
    C = 3.0 * (gamma + 1.0) / 16.0
    
    # 5. Print the final equation with the calculated numerical coefficients.
    print("The analytical solution for the normalized density profile is:")
    print(f"rho(x)/rho_0 = {A:.2f} + ({B:.2f}) * tanh({C:.2f} * x/L)")
    print("\nSimplified form:")
    # Using abs(B) to present it in the standard form A - B*tanh(...)
    print(f"rho(x)/rho_0 = {A:.2f} - {abs(B):.2f} * tanh({C:.2f} * x/L)")

# Execute the function to get the solution
solve_density_profile()
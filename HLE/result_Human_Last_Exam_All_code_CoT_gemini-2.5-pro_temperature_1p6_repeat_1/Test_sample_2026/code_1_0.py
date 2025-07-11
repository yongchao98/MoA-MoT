import math

def solve_problem():
    """
    This function calculates the area of the surface defined by the initial conditions
    (x0, y0, z0) for which the given nonlinear boundary value problem has solutions.
    """
    
    # As derived from the solvability conditions, the relationship between
    # the initial values (x0, y0, z0) and the parameters (alpha, beta, gamma) is:
    #   alpha = (y0^2 + z0^2) * (1 - e^-T)
    #   beta  = (x0^2 + z0^2) * (1 - e^-T)
    #   gamma = (x0^2 + y0^2) * (1 - e^-T)
    
    # Summing these three equations gives:
    #   alpha + beta + gamma = 2 * (x0^2 + y0^2 + z0^2) * (1 - e^-T)
    
    # The problem provides the constraint:
    #   alpha + beta + gamma = 10^25 * (1 - e^-T)
    alpha_beta_gamma_sum_coeff = 10**25
    
    # By substituting the constraint into the summed equation, the term (1 - e^-T) cancels out:
    #   10^25 * (1 - e^-T) = 2 * (x0^2 + y0^2 + z0^2) * (1 - e^-T)
    #   10^25 = 2 * (x0^2 + y0^2 + z0^2)
    
    # This simplifies to the equation of a sphere:
    #   x0^2 + y0^2 + z0^2 = 10^25 / 2
    
    # The squared radius (R^2) of this sphere is:
    R_squared_val = alpha_beta_gamma_sum_coeff / 2
    
    # The surface area of a sphere is given by the formula A = 4 * pi * R^2.
    area = 4 * math.pi * R_squared_val
    
    # --- Output the derivation and result ---
    print("The set of initial values (x0, y0, z0) for which solutions exist forms a sphere.")
    print("The derivation leads to the following equation for this sphere:")
    
    # Print the equation with numbers
    # The final equation is x0^2 + y0^2 + z0^2 = R^2
    print(f"\nx0^2 + y0^2 + z0^2 = {alpha_beta_gamma_sum_coeff} / {2}")
    
    # The equation for the surface area is A = 4 * pi * R^2
    print(f"\nThe area of this sphere is A = 4 * pi * ({alpha_beta_gamma_sum_coeff} / {2})")
    
    # Simplify the equation for the area
    # A = 2 * pi * 10^25
    print(f"This simplifies to the final equation for the area: A = {2} * pi * {alpha_beta_gamma_sum_coeff}")
    
    print(f"\nCalculated numeric value for the area:")
    print(f"{area:.4e}")

solve_problem()
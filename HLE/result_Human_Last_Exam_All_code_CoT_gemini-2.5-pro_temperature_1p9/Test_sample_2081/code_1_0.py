import math

def solve_for_R():
    """
    Calculates the radius R of the sphere of initial values for which solutions exist.
    """
    # Given T, we can calculate its numerical value
    # T = ln(10^34)
    T = math.log(10**34)
    
    # Calculate e^T and e^(2T)
    # Due to floating point precision, e^T will be exactly 10^34
    e_T = 10**34
    e_2T = e_T * e_T # To avoid potential overflow issues with math.exp(2*T) on some platforms
    
    # Calculate alpha using the given formula
    alpha = 0.5 * (e_2T - 1)
    
    # The set of initial values (x_0, y_0, z_0) is defined by the solvability condition,
    # which describes a sphere:
    # 1 * x_0^2 + 1 * y_0^2 + 1 * z_0^2 = R^2
    
    # We found that R^2 = alpha / (1 - e^(-T))
    
    # Let's calculate the components of the R^2 formula
    one_minus_exp_neg_T = 1 - math.exp(-T)
    
    # Now, calculate R^2
    R_squared = alpha / one_minus_exp_neg_T

    # The final equation for the initial values is a sphere.
    # We output the coefficients and the constant term (R^2).
    coeff_x_sq = 1
    coeff_y_sq = 1
    coeff_z_sq = 1

    print("The final equation for the initial values (x_0, y_0, z_0) is that of a sphere:")
    print(f"{coeff_x_sq}*x_0^2 + {coeff_y_sq}*y_0^2 + {coeff_z_sq}*z_0^2 = R^2")
    print("\nCalculated numbers for the equation:")
    print(f"Coefficient for x_0^2: {coeff_x_sq}")
    print(f"Coefficient for y_0^2: {coeff_y_sq}")
    print(f"Coefficient for z_0^2: {coeff_z_sq}")
    print(f"Value of R^2: {R_squared}")
    
    # Finally, calculate R by taking the square root
    R = math.sqrt(R_squared)
    
    print(f"\nThe radius of the sphere is R = {R}")

solve_for_R()
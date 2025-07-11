import math

def solve_problem():
    """
    This function solves the problem by calculating the radius R based on the solvability condition
    of the given nonlinear boundary-value problem.
    """
    print("Step 1: Deriving the equation for R.")
    print("------------------------------------")
    print("The solvability condition for the given weakly nonlinear problem leads to a constraint on the initial values (x0, y0, z0) of the unperturbed solution.")
    print("This condition is: (x0^2 + y0^2 + z0^2) * (1 - e^(-T)) = alpha")
    print("This equation describes a sphere in the space of initial values.")
    print("The radius R of this sphere is what we need to find.")
    print("So, R^2 is given by the equation: R^2 = alpha / (1 - e^(-T))")
    print("\n")

    print("Step 2: Substituting the given values for T and alpha.")
    print("-----------------------------------------------------")
    # T = ln(10^34), so e^T = 10^34 and e^(-T) = 10^(-34)
    # alpha = (1/2) * (e^(2T) - 1) = 0.5 * ((10^34)^2 - 1)
    e_T = 10**34
    e_neg_T = 10**-34
    alpha_numerator = e_T**2 - 1
    alpha_denominator = 2

    denominator = 1 - e_neg_T
    
    print(f"Given T = ln(10^34), we have e^T = {e_T} and e^(-T) = {e_neg_T:.1e}.")
    print(f"Given alpha = (1/2) * (e^(2T) - 1), we have alpha = ({alpha_numerator}) / {alpha_denominator}.")
    
    print("\nPlugging these into the equation for R^2:")
    print(f"R^2 = (({alpha_numerator}) / {alpha_denominator}) / (1 - {e_neg_T:.1e})")
    print("\n")

    print("Step 3: Simplifying the equation and calculating R.")
    print("-------------------------------------------------")
    # The expression for R^2 can be simplified algebraically:
    # R^2 = (0.5 * (e_T^2 - 1)) / (1 - 1/e_T)
    #     = (0.5 * (e_T - 1) * (e_T + 1)) / ((e_T - 1) / e_T)
    #     = 0.5 * (e_T + 1) * e_T
    # This simplified form is better for numerical computation.
    
    # We use Python's arbitrary-precision integers for an exact calculation of R^2.
    R_squared_numerator_term1 = e_T + 1
    R_squared_numerator_term2 = e_T
    R_squared_denominator = 2
    
    print("The simplified final equation for R^2 is:")
    # We display the numbers in the equation as requested.
    print(f"R^2 = (({R_squared_numerator_term1}) * ({R_squared_numerator_term2})) / {R_squared_denominator}")

    R_squared_val = (R_squared_numerator_term1 * R_squared_numerator_term2) // R_squared_denominator
    R_val = math.sqrt(R_squared_val)

    print("\nFinal calculation:")
    print(f"R = sqrt({R_squared_val})")
    print(f"R = {R_val:.15e}")

solve_problem()
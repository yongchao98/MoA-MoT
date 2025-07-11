import math

def solve_for_t0():
    """
    Calculates the positive value of t_0 based on the solvability condition
    of the given boundary-value problem.
    """
    # Given parameters
    alpha = 10**16
    R = math.log(100.0 / 99.0)

    # The solvability condition leads to the algebraic equation:
    # t_0^2 * Integral_1 = Integral_2
    # where:
    # Integral_1 = integral from 0 to R of integral from 0 to inf of (phi^2 * v) dx dt
    # Integral_2 = integral from 0 to inf of (alpha * v(x,0)) dx
    #
    # With phi = exp(-x+t) and v = exp(-x-t), we have:
    # Integral_1 = (exp(R) - 1) / 3
    # Integral_2 = alpha
    #
    # So, the equation is: t_0^2 * ( (exp(R) - 1) / 3 ) = alpha

    # Calculate the terms of the equation
    exp_R = math.exp(R)
    integral_1_val = (exp_R - 1) / 3
    integral_2_val = alpha

    # The final equation for t_0^2
    # t_0^2 = 3 * alpha / (exp(R) - 1)
    t0_squared = 3 * alpha / (exp_R - 1)

    # t_0 must be positive
    t0 = math.sqrt(t0_squared)

    # Print the equation with the numerical values
    print("The derived algebraic equation for t_0 is:")
    print(f"(t_0)^2 * ( (exp({R}) - 1) / 3 ) = {alpha}")
    print("Substituting the calculated values:")
    print(f"(t_0)^2 * {integral_1_val} = {integral_2_val}")
    print("\nSolving for t_0:")
    print(f"t_0 = sqrt({integral_2_val} / {integral_1_val})")
    print(f"t_0 = sqrt({t0_squared})")
    
    # Print the final answer
    print("\nThe positive value of t_0 is:")
    print(t0)
    
    # Return the final answer in the required format
    print(f"\n<<<{t0}>>>")

solve_for_t0()
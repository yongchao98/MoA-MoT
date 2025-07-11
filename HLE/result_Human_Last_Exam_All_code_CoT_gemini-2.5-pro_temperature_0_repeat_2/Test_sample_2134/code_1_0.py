import math

def solve_trajectory():
    """
    Calculates the particle's position x(t) at t = 2*sqrt(3).
    """
    # Based on the analysis, the particle's trajectory is given by x(t) = x(0) + t^2 / 4.

    # Step 1: Define the initial position x(0).
    # The expression is x(0) = 3 + (6*(3-sqrt(3)))^(1/3) + (6*(3+sqrt(3)))^(1/3)
    # This is equivalent to x(0) = 3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)

    sqrt_3 = math.sqrt(3)
    
    # The terms inside the cube roots
    term_in_cbrt_1 = 18 - 6 * sqrt_3
    term_in_cbrt_2 = 18 + 6 * sqrt_3

    # The two radical terms in the expression for x(0)
    radical_term_1 = term_in_cbrt_1**(1/3)
    radical_term_2 = term_in_cbrt_2**(1/3)

    # The constant part of x(0)
    constant_part_x0 = 3
    
    # Calculate x(0)
    x_0 = constant_part_x0 + radical_term_1 + radical_term_2

    # Step 2: Calculate the displacement term t^2 / 4.
    t = 2 * sqrt_3
    displacement = t**2 / 4

    # Step 3: Calculate the final position x(t).
    x_t = x_0 + displacement

    # Output the numbers in the final equation as requested.
    # The final equation is x(t) = x(0) + 3, which expands to:
    # x(t) = (3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)) + 3
    # x(t) = 6 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)
    
    final_constant_part = constant_part_x0 + displacement

    print("The final position x(t) is given by the equation:")
    print(f"x(t) = {final_constant_part} + ({term_in_cbrt_1})^(1/3) + ({term_in_cbrt_2})^(1/3)")
    print("\nBreaking it down into numerical values:")
    print(f"x(t) = {final_constant_part} + {radical_term_1} + {radical_term_2}")
    print(f"x(t) = {x_t}")

solve_trajectory()
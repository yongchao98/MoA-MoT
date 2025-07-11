import math

def solve_trajectory():
    """
    Calculates the particle's position at t = 2*sqrt(3)
    based on the provided Bohmian mechanics problem.
    """
    # Step 1: Calculate the numerical value of the initial position x(0).
    # x(0) = 3 + (6*(3 - sqrt(3)))^(1/3) + (6*(3 + sqrt(3)))^(1/3)
    sqrt_3 = math.sqrt(3)
    
    term_inside_cbrt1 = 6 * (3 - sqrt_3)
    term_inside_cbrt2 = 6 * (3 + sqrt_3)
    
    cbrt_term1 = term_inside_cbrt1**(1/3)
    cbrt_term2 = term_inside_cbrt2**(1/3)
    
    x0 = 3 + cbrt_term1 + cbrt_term2
    
    # Step 2: The trajectory of the particle is x(t) = t^2/2 + x(0).
    # We need to find the position at t = 2*sqrt(3).
    t = 2 * sqrt_3
    
    # The time-dependent displacement is t^2/2.
    time_displacement = t**2 / 2
    
    # Step 3: Calculate the final position.
    final_position = time_displacement + x0
    
    # Step 4: Print the final equation with each numerical component.
    # The final equation is x(2*sqrt(3)) = 6 + x(0)
    print("The final position is the sum of the time-dependent displacement and the initial position.")
    print(f"Final Position = (Time Displacement) + (Initial Position)")
    print(f"{final_position:.4f} = {time_displacement:.4f} + {x0:.4f}")

solve_trajectory()
import math

def solve_trajectory():
    """
    Calculates the final position of the particle based on the derived trajectory equation.
    """
    # The final equation for the particle's position is x(t) = x(0) - 3.
    # We need to calculate the initial position x(0) and the final position x(t).

    # Part 1: Calculate the value of x(0).
    # x(0) = 3 + ∛(6(3-√3)) + ∛(6(3+√3))
    # Let's calculate the complex part of x(0).
    y = (6 * (3 - math.sqrt(3)))**(1/3) + \
        (6 * (3 + math.sqrt(3)))**(1/3)
    
    # The full initial position x(0)
    x0 = 3 + y
    
    # Part 2: Calculate the displacement at t = 2√3.
    # The displacement is -t^2/4.
    t = 2 * math.sqrt(3)
    displacement = -t**2 / 4

    # Part 3: Calculate the final position x(t).
    # x(t) = x(0) + displacement
    x_final = x0 + displacement

    # Output the numbers in the final equation as requested.
    print(f"Initial position x(0): {x0}")
    print(f"Displacement at t=2*sqrt(3): {displacement}")
    print(f"Final position x(t) = x(0) + displacement = {x0} + ({displacement})")
    print(f"The calculated final position is: {x_final}")

# Run the solver
solve_trajectory()
import math

def solve_particle_emitter_problem():
    """
    Calculates the minimum ratio of the cube of the surface area to the
    square of the volume of the region traversed by the particles.

    The analytical solution shows the minimum ratio is a constant value,
    given by the expression: 36 * pi * (3 + 2 * sqrt(3)).
    """

    # The components of the final analytical solution
    factor_1 = 36
    factor_2 = 3
    factor_3 = 2
    sqrt_operand = 3

    pi_val = math.pi
    sqrt_3_val = math.sqrt(sqrt_operand)

    # Print the equation that gives the minimum ratio
    # This fulfills the requirement to output each number in the final equation.
    print("The final equation for the minimum ratio (Q) is:")
    print(f"Q = {factor_1} * pi * ({factor_2} + {factor_3} * sqrt({sqrt_operand}))")
    print("\nCalculating the numerical value...")

    # Calculate the numerical result
    final_ratio = factor_1 * pi_val * (factor_2 + factor_3 * sqrt_3_val)

    # Print the result
    print(f"\nThe minimum ratio is: {final_ratio}")

# Execute the function to solve the problem
solve_particle_emitter_problem()
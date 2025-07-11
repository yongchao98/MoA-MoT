import math

def solve_gravity_shape_problem():
    """
    Calculates the maximum distance from point A for a shape of uniform density
    and 1 cubic meter volume that creates the strongest gravitational field at A.

    The shape is given by r(θ) = C * sqrt(cos(θ)).
    The furthest point is at θ=0, where the distance is C.
    The volume constraint V=1m³ leads to C = (15 / (4 * π))^(1/3).
    """

    # Define the constants from the derived formula
    numerator = 15.0
    denominator_factor = 4.0
    exponent = 1.0 / 3.0

    print("The problem is to find the furthest point on a 1m³ shape that maximizes the gravitational field at a point A.")
    print("The distance of this furthest point is given by the formula: C = (numerator / (denominator_factor * π))^(exponent)")
    print(f"Here, numerator = {numerator}, denominator_factor = {denominator_factor}, and exponent = 1/3.")
    print("-" * 50)

    # --- Step-by-step calculation ---

    # Step 1: Calculate the value of pi
    pi_val = math.pi
    print(f"Step 1: Using the value of π ≈ {pi_val}")

    # Step 2: Calculate the denominator
    denominator = denominator_factor * pi_val
    print(f"\nStep 2: Calculate the denominator ({denominator_factor} * π)")
    print(f"   {denominator_factor} * {pi_val} = {denominator}")

    # Step 3: Calculate the base of the power
    base = numerator / denominator
    print(f"\nStep 3: Calculate the term inside the parenthesis ({numerator} / {denominator})")
    print(f"   {numerator} / {denominator} = {base}")

    # Step 4: Calculate the final result
    result = base ** exponent
    print(f"\nStep 4: Calculate the cube root of the result ({base}^(1/3))")
    print(f"   {base} ** (1/3) = {result}")

    print("-" * 50)

    # --- Final Equation Output ---
    print("The final equation with all numbers shown is:")
    print(f"({numerator} / ({denominator_factor} * {pi_val})) ** (1/3) = {result}")
    print(f"\nThe furthest point on the surface of the playdough is {result} meters from point A.")


if __name__ == '__main__':
    solve_gravity_shape_problem()
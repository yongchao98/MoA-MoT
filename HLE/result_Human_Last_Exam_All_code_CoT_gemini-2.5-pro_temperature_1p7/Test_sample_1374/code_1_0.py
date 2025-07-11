import math

def solve_gravity_problem():
    """
    Calculates the furthest point on a 1 m^3 volume of uniform density,
    shaped to produce the strongest possible gravitational field at a point A.

    The shape that maximizes the gravitational field is bounded by the surface
    r^2 = C * cos(theta).
    The volume of this shape is V = (4 * pi * C^(3/2)) / (3 * (5/2)) = (8 * pi * C^(3/2)) / 15.
    Wait, the derivation in thought was V = 4*pi / (15*k^(3/2)) with k=1/C. Let's recheck.
    V = ∫[0, 2π] dφ ∫[0, π/2] sin(θ) dθ ∫[0, sqrt(C*cos(θ))] r^2 dr
      = 2π ∫[0, π/2] sin(θ) * [r^3/3]_0^sqrt(C*cos(θ)) dθ
      = (2π/3) * C^(3/2) ∫[0, π/2] cos^(3/2)(θ)sin(θ) dθ
    Let u = cos(θ), du = -sin(θ)dθ. Integral becomes ∫[1, 0] u^(3/2) (-du) = ∫[0, 1] u^(3/2) du = [u^(5/2)/(5/2)]_0^1 = 2/5.
    So, V = (2π/3) * C^(3/2) * (2/5) = (4π/15) * C^(3/2). This is correct.
    
    Given V = 1 m^3:
    1 = (4 * pi / 15) * C^(3/2) => C^(3/2) = 15 / (4 * pi) => C = (15 / (4 * pi))^(2/3).
    
    The furthest distance r_max occurs at theta = 0, where r^2 = C * cos(0) = C.
    So, r_max = sqrt(C) = sqrt((15 / (4 * pi))^(2/3)) = (15 / (4 * pi))^(1/3).
    This confirms the final formula.
    """
    
    # Constants from the final derived equation for r_max
    numerator = 15
    denominator_factor = 4
    pi_val = math.pi
    exponent_numerator = 1
    exponent_denominator = 3

    # The equation to be solved
    print("The problem is to find the maximum distance from a point A to the surface of a shape")
    print("that has a volume of 1 cubic meter and maximizes the gravitational field at A.")
    print("This distance (r_max) can be calculated with the following equation:")
    print(f"r_max = ({numerator} / ({denominator_factor} * pi)) ^ ({exponent_numerator}/{exponent_denominator})\n")

    # Outputting each number in the final equation
    print("Here are the numbers used in the equation:")
    print(f"Numerator: {numerator}")
    print(f"Denominator factor: {denominator_factor}")
    print(f"Value of pi: {pi_val}")
    print(f"Exponent: {exponent_numerator}/{exponent_denominator}\n")

    # Calculate the result
    base = numerator / (denominator_factor * pi_val)
    r_max = base ** (exponent_numerator / exponent_denominator)

    # Print the final result
    print(f"The furthest distance from A is: {r_max} meters.")

solve_gravity_problem()
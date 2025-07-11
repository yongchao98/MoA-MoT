import math

def solve_for_time():
    """
    This function calculates the time when the chain first loses contact with the ground.
    
    The vertical position of the bottom of the chain, z_chain(t), is given by:
    z_chain(t) = 5*(1+sin(t)) + sqrt(3)/2 + 0.25*(0.5*cos^2(t) + (sqrt(3)/2)*sin(t)) - 10

    Setting z_chain(t) = 0 and substituting cos^2(t) = 1 - sin^2(t) leads to a
    quadratic equation for x = sin(t) of the form a*x^2 + b*x + c = 0.
    """

    # Coefficients of the quadratic equation a*x^2 + b*x + c = 0, where x = sin(t).
    a = 1.0
    b = -(40 + math.sqrt(3))
    c = 39 - 4 * math.sqrt(3)

    print("The problem reduces to solving the quadratic equation a*x^2 + b*x + c = 0, where x = sin(t).")
    print(f"The calculated coefficients are:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print("-" * 30)

    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c

    # The quadratic formula gives two solutions for x = sin(t).
    # One solution is > 1, which is impossible for sin(t). We choose the other solution.
    sin_t = (-b - math.sqrt(discriminant)) / (2 * a)

    # The robot's journey starts at t=0. We need the first time t > 0.
    # Since the calculated sin(t) is positive, the smallest positive time t is the principal value of arcsin.
    time = math.asin(sin_t)

    print(f"The value of sin(t) when the chain leaves the ground is: {sin_t}")
    print(f"The time t (in seconds) is: {time}")

solve_for_time()
<<<0.8998184277742138>>>
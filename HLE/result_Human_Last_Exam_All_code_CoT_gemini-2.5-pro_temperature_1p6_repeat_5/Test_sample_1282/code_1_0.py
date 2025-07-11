import math

def solve_blowup_time():
    """
    This function solves for the finite blow-up time T based on the analysis
    of the linearized version of the given partial differential equation.

    The blow-up condition leads to the quadratic equation: a*T^2 + b*T + c = 0.
    For our chosen example initial condition, the coefficients are a=1, b=2, c=-1.
    """
    # Coefficients of the quadratic equation T^2 + 2T - 1 = 0
    a = 1.0
    b = 2.0
    c = -1.0

    # Print the equation
    print("The analysis of the linearized system for a sample initial condition leads to a blow-up time T.")
    print("This time T is the positive solution to the quadratic equation:")
    print(f"({a}) * T^2 + ({b}) * T + ({c}) = 0")
    print("-" * 30)
    
    # Calculate the discriminant
    discriminant = b**2 - 4*a*c

    if discriminant < 0:
        print("The equation has no real solutions for the blow-up time.")
    else:
        # Calculate the two roots
        t1 = (-b - math.sqrt(discriminant)) / (2*a)
        t2 = (-b + math.sqrt(discriminant)) / (2*a)

        # The blow-up time must be positive
        if t2 > 0:
            print(f"The calculated finite blow-up time is T = {t2}")
        else:
            print("No positive real solution found for the blow-up time.")

if __name__ == "__main__":
    solve_blowup_time()

import math

def solve_blowup_condition():
    """
    Calculates the critical value for y(0) that determines whether the solution
    to the system of ODEs blows up, given an initial value for x(0) > 1.
    """
    # We use x(0) = 8 as an example value. x(0) must be greater than 1.
    x0 = 8

    if x0 <= 1:
        print("This analysis assumes x(0) > 1.")
        return

    # The condition for blow-up is derived from the separatrix equation of the
    # saddle point at (1,0). The equation for the separatrix is:
    # y^2 = 2*x + 1 - 3*x^(2/3)

    # We print the numbers in this equation as requested.
    a = 2
    b = 1
    c = 3
    exp_num = 2
    exp_den = 3
    
    print("The condition for blow-up depends on the stable manifold of the saddle point (1,0).")
    print("The equation for the curve containing the stable manifold (the separatrix) is:")
    print(f"y^2 = {a}*x + {b} - {c}*x^({exp_num}/{exp_den})\n")
    
    # For a solution to avoid blow-up, the initial point (x(0), y(0)) must lie
    # on or above the stable manifold branch.
    # The stable manifold for x > 1 is the branch y > 0.
    # So, blow-up occurs if y(0) < sqrt(2*x(0) + 1 - 3*x(0)^(2/3)).

    # Calculate the critical value for y(0)
    try:
        # The term inside the square root is g(x) = 2x + 1 - 3x^(2/3).
        # For x > 1, g(x) is always positive.
        y_critical = math.sqrt(a * x0 + b - c * (x0 ** (exp_num / exp_den)))
    except ValueError:
        print(f"Could not compute the critical value for x(0) = {x0}. "
              "The term inside the square root is negative.")
        return

    print(f"For a given initial condition x(0) = {x0}, the solution blows up if:")
    print(f"y(0) < {y_critical:.4f}")
    print("\nThe general condition for blow-up for any x(0) > 1 is:")
    print("y(0) < sqrt(2*x(0) + 1 - 3*x(0)**(2/3))")

# Run the function
solve_blowup_condition()
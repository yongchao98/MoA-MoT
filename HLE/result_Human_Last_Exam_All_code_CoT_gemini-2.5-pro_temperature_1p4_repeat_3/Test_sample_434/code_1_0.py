import math

def find_blowup_condition(x0):
    """
    For the given system of ODEs, this function calculates the critical value for y(0) 
    based on an initial value x(0) > 1, determining the range for which the solution blows up.

    Args:
        x0 (float): The initial condition for x(t), which must be greater than 1.
    """
    if x0 <= 1:
        print("Error: The analysis is only valid for x(0) > 1.")
        print("Please provide a value for x(0) that is greater than 1.")
        return

    # These are the numeric coefficients and powers from the derived formula for the blow-up condition.
    c1 = 1
    c2 = 2
    c3 = 3
    p_num = 2
    p_den = 3

    # The condition for blow-up is y(0) < sqrt(1 + 2*x(0) - 3*x(0)^(2/3)).
    # We first calculate the expression inside the square root.
    critical_y_squared = c1 + c2 * x0 - c3 * x0**(p_num / p_den)
    
    # As derived, for x0 > 1, critical_y_squared is always positive.
    critical_y_bound = math.sqrt(critical_y_squared)

    # First, print the general inequality showing all the numbers involved.
    print("For the solution to blow up, the initial condition y(0) must satisfy the inequality:")
    print(f"y(0) < sqrt({c1} + {c2}*x(0) - {c3}*x(0)**({p_num}/{p_den}))")
    
    # Then, print the result for the specific value of x(0) provided.
    print(f"\nFor the specific case where x(0) = {x0}, the condition is:")
    print(f"y(0) < {critical_y_bound:.6f}")

# We are given that x(0) > 1. Let's use x(0) = 4 as an example.
x_initial = 4.0
find_blowup_condition(x_initial)
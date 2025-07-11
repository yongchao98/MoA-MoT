import math

def blowup_condition(x0):
    """
    Calculates the threshold for y(0) for a given x(0) > 1 that leads to a blow-up.
    
    Args:
        x0 (float): The initial condition x(0), must be greater than 1.
        
    Returns:
        None. Prints the condition for y(0).
    """
    if x0 <= 1:
        print("Error: x(0) must be greater than 1.")
        return

    # The constants in the threshold equation y(0) <= sqrt(c2*x(0) + c1 - c3*x(0)**(p1/p2))
    c2, c1, c3 = 2, 1, 3
    p1, p2 = 2, 3
    
    # Calculate the term inside the square root
    inner_term = c2 * x0 + c1 - c3 * (x0**(p1/p2))
    
    # The square root is well-defined for x0 > 1
    threshold = math.sqrt(inner_term)
    
    print("For the solution to blow up, the initial condition y(0) must satisfy the inequality:")
    print(f"y(0) <= sqrt({c2}*x(0) + {c1} - {c3}*x(0)^({p1}/{p2}))")
    print("\nFor the specific case where x(0) =", x0, ", the condition is:")
    print(f"y(0) <= {threshold:.4f}")

# Example usage with x(0) = 2
x_initial = 2
blowup_condition(x_initial)

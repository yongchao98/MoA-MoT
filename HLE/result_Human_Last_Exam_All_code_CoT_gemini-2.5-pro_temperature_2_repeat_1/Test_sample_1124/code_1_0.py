import math

def solve_modified_logistic_map():
    """
    This function modifies the standard logistic map and finds its new equilibrium point.
    The modification adds the term (R - X) / 2 to the standard map.
    """
    R = 3.57

    # The modified map is X_n+1 = R * X_n * (1 - X_n) + (R - X_n) / 2.
    # An equilibrium point X* satisfies X* = R*X*(1-X*) + (R-X*)/2.
    # Rearranging this gives the quadratic equation:
    # 2*R*X**2 + (3 - 2*R)*X - R = 0
    a = 2 * R
    b = 3 - 2 * R
    c = -R

    # Solve the quadratic equation using the quadratic formula
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        print("The equilibrium point is not a real number.")
        return
        
    # We are interested in the positive root, as it is close to 1.05
    X_star = (-b + math.sqrt(discriminant)) / (2 * a)

    # Now we verify that this X_star is indeed the equilibrium point.
    # We calculate the left-hand side (LHS) and right-hand side (RHS) of the equation.
    lhs = X_star
    rhs = R * X_star * (1 - X_star) + (R - X_star) / 2

    print("The modified logistic map is: X_n+1 = R * X_n * (1 - X_n) + (R - X_n) / 2")
    print(f"\nFor the parameter R = {R}, the new equilibrium point is calculated to be approximately {X_star:.4f}.")
    print("\nTo verify this, we check if the point satisfies the equilibrium condition X = g(X, R).")
    print("Plugging the values into the equation:")
    # Print the equation with all the numbers, as requested.
    print(f"EQUATION: {X_star:.6f} = {R} * {X_star:.6f} * (1 - {X_star:.6f}) + ({R} - {X_star:.6f}) / 2")

    # Print the evaluated sides of the equation to show they are equal.
    print(f"\nCalculated LHS: {lhs:.6f}")
    print(f"Calculated RHS: {rhs:.6f}")
    
solve_modified_logistic_map()
# The value for the answer box should be the equilibrium point
a = 2 * 3.57
b = 3 - 2 * 3.57
c = -3.57
discriminant = b**2 - 4*a*c
X_star = (-b + math.sqrt(discriminant)) / (2 * a)
print(f'<<<{X_star:.2f}>>>')

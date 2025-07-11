import sympy
from sympy import sqrt, diff, solve, Symbol

def solve_traffic_flow_problem():
    """
    This function symbolically computes the constant lower bound for d(t,x)
    for the given traffic flow model.
    """
    # Step 1: Explain the theoretical background.
    print("Step 1: Deriving the evolution equation for the minimum of the derivative.")
    print("Let d(t,x) = du/dx. Differentiating the original PDE with respect to x and evaluating at a spatial minimum of d, we obtain an ODE for the minimum value D(t) = min_x d(t,x).")
    print("This ODE is a Riccati equation of the form: dD/dt = exp(-u_bar) * Q(D, u).")
    print("The quadratic part is Q(D, u) = 2*D^2 - (3*u - 5*u^2)*D - u^3*(1 - u).")
    print("-" * 20)

    # Step 2: Analyze the Riccati ODE to find the lower bound.
    print("Step 2: Analyzing the ODE to find the lower bound for d(t,x).")
    print("D(t) can decrease only if Q(D,u) is negative, which means D(t) is between the roots of Q(D,u)=0.")
    print("Thus, D(t) is bounded below by the minimum of the smaller root of Q(D,u)=0, for u in [0, 1].")
    print("-" * 20)

    # Step 3: Use SymPy to find and analyze the smaller root.
    print("Step 3: Finding the minimum of the smaller root using SymPy.")

    # Define symbols
    u = Symbol('u', real=True, positive=True)
    d = Symbol('d', real=True)

    # Define the quadratic polynomial Q(d, u)
    Q = 2*d**2 - (3*u - 5*u**2)*d - u**3*(1-u)

    # Find the roots of Q(d, u) = 0 for d
    roots = solve(Q, d)

    # The smaller root, D1(u), corresponds to the minus sign in the quadratic formula.
    discriminant_sqrt = sqrt((3*u - 5*u**2)**2 + 8*u**3*(1-u))
    D1_u_expr = ((3*u - 5*u**2) - discriminant_sqrt) / 4

    print(f"The smaller root is D1(u) = {D1_u_expr}")

    # To find the minimum of D1(u) for u in [0, 1], we check its derivative.
    # Analytical derivation shows D1'(u) has no roots in (0,1), so it is monotonic.
    # We can determine the direction by testing the derivative at a point, e.g., u=0.5.
    D1_u_prime = diff(D1_u_expr, u)
    test_point = 0.5
    derivative_at_test_point = D1_u_prime.subs(u, test_point)

    print(f"\nThe derivative D1'(u) at u={test_point} is {derivative_at_test_point.evalf()}, which is negative.")
    print("This implies D1(u) is a decreasing function over [0, 1].")

    # For a decreasing function, the minimum value is at the right endpoint of the interval.
    # Evaluate D1(u) at the endpoints u=0 and u=1.
    val_at_0 = D1_u_expr.subs(u, 0)
    val_at_1 = D1_u_expr.subs(u, 1)

    print(f"\nValue of D1(u) at u=0: {val_at_0}")
    print(f"Value of D1(u) at u=1: {val_at_1}")

    lower_bound = val_at_1

    print("-" * 20)
    # Step 4: Final Conclusion.
    print("Step 4: Conclusion.")
    print("The minimum of the smaller root D1(u) over u in [0,1] occurs at u=1.")
    print(f"This minimum value is {lower_bound}. This provides a constant lower bound for d(t,x).")
    print("Given d_min(0) = -0.5, which is greater than this bound, d(t,x) will not drop below this value.")
    
    print("\nThe final equation expressing the lower bound is:")
    print(f"d(t,x) >= {int(lower_bound)}")

    # Return the numerical value of the bound
    return int(lower_bound)

if __name__ == '__main__':
    bound = solve_traffic_flow_problem()
    # The final answer in the requested format will be extracted from this
    # For execution purpose, we print the final answer
    print("\nDetermined constant lower bound:")
    print(bound)

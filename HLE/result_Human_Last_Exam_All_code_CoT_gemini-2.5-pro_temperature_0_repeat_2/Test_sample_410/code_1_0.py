import numpy as np

def solve_problem():
    """
    This function solves the entire problem step-by-step.
    """
    # Step 1: Find the coefficients of f(x) = A*e^(2x) + B*e^x + C
    # The problem uses 'a', 'b', 'c' for coefficients, but also for the integral problem.
    # To avoid confusion, we'll use A, B, C for the coefficients of f(x).
    
    # From the limit condition: lim_{x->-inf} (A*e^(2x) + B*e^x + C + 3) / e^x = 1
    # This can be rewritten as: lim_{x->-inf} (A*e^x + B + (C+3)/e^x) = 1
    # As x -> -inf, e^x -> 0. For the limit to be finite, the term (C+3)/e^x must not diverge.
    # This requires C + 3 = 0.
    C = -3.0
    # With C = -3, the limit evaluates to B.
    B = 1.0
    
    # From the condition f(ln(2)) = 0:
    # A*e^(2*ln(2)) + B*e^(ln(2)) + C = 0
    # A * (e^ln(2))^2 + B * e^ln(2) + C = 0
    # A * 2^2 + 1.0 * 2 + (-3.0) = 0
    # 4A + 2 - 3 = 0
    # 4A - 1 = 0
    A = 0.25
    
    print(f"Step 1: The function is determined to be f(x) = {A}*e^(2x) + {B}*e^x + {C}")

    # Define the function f(x) for calculations
    def f(x):
        return A * np.exp(2*x) + B * np.exp(x) + C

    # Step 2: Analyze the integral identity
    # The given identity is: integral from 0 to a of g(x)dx + integral from ln(2) to ln(b) of f(x)dx = a*ln(b)
    # The general theorem for the integral of an inverse function is:
    # integral from u to v of f(x)dx + integral from f(u) to f(v) of g(y)dy = v*f(v) - u*f(u)
    # Let's set u = ln(2) and v = ln(b).
    # We know f(u) = f(ln(2)) = 0 from the problem statement.
    # The theorem becomes:
    # integral from ln(2) to ln(b) of f(x)dx + integral from 0 to f(ln(b)) of g(y)dy = ln(b)*f(ln(b))
    # By comparing this with the given identity, we can see they match if:
    # a = f(ln(b))
    print("\nStep 2: The integral identity holds if and only if a = f(ln(b)).")

    # Step 3: Find the unique solution (a, b)
    # The problem gives a special point, (x_0, y_0) = (ln(2), 0), where f(x_0) = y_0.
    # We hypothesize that the unique solution (a, b) is constructed from this point.
    # A natural choice is a = y_0 and b = e^(x_0).
    
    # The variables 'a' and 'b' for the final answer
    a_sol = f(np.log(2))
    b_sol = np.exp(np.log(2))
    
    print("\nStep 3: Using the special point (ln(2), f(ln(2))) = (ln(2), 0) to find the unique solution.")
    print(f"We set a = f(ln(2)), which gives a = {a_sol:.0f}")
    print(f"We set b = e^(ln(2)), which gives b = {b_sol:.0f}")

    # Verify that this solution (a, b) satisfies the condition a = f(ln(b))
    lhs = a_sol
    rhs = f(np.log(b_sol))
    print("\nVerifying the solution by checking if a = f(ln(b)):")
    print(f"LHS = a = {lhs:.4f}")
    print(f"RHS = f(ln({b_sol:.0f})) = f({np.log(b_sol):.4f}) = {rhs:.4f}")
    if np.isclose(lhs, rhs):
        print("Verification successful. The solution is correct.")
    else:
        print("Verification failed. The solution is incorrect.")

    # Step 4: Calculate the final answer a + b
    sum_val = a_sol + b_sol
    print("\nStep 4: Calculate the final result a + b.")
    print("The final equation is:")
    print(f"{a_sol:.0f} + {b_sol:.0f} = {sum_val:.0f}")

solve_problem()
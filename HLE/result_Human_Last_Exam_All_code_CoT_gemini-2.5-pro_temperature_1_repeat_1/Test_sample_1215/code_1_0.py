def solve_logic_puzzle():
    """
    Solves the propositional logic puzzle by determining the number of essential variables.
    """
    
    # The problem asks for the minimum number of variables in a formula psi,
    # which is equivalent to a formula phi with n variables.
    # This minimum number is the number of 'essential' variables in phi.
    # Let's denote the number of essential variables by k.

    # Let n be the number of atomic variables in phi (n >= 2).
    # Let k be the number of essential variables in phi (1 <= k <= n).
    # Let m be the number of dummy variables, so m = n - k.
    
    # Condition 1: phi has 2^(n-1) satisfying assignments (models).
    # The number of models can also be expressed as:
    # (Number of satisfying assignments for essential vars) * 2^(number of dummy vars)
    # Let N_k_true be the number of satisfying assignments for the k essential variables.
    # The total number of models is N_k_true * 2^m = N_k_true * 2^(n-k).
    
    # So, we have the equation:
    # N_k_true * 2^(n-k) = 2^(n-1)
    # Solving for N_k_true:
    # N_k_true = 2^(n-1) / 2^(n-k) = 2^((n-1) - (n-k)) = 2^(k-1).
    
    # This equation holds for any k from 1 to n. For example:
    # If k=1, N_1_true = 2^(1-1) = 1. A function of 1 variable with 1 model is f(p)=p. This is possible.
    # If k=n, N_n_true = 2^(n-1). A function of n variables with 2^(n-1) models is the parity function (p1 XOR p2 XOR ...). This is also possible.
    
    # The key to solving the problem lies in the word "discernible" from Condition 1:
    # "It has exactly 2^(n-1) discernible truth-value assignments."
    # This non-standard term implies a special property. A logical interpretation is that
    # any two distinct models of phi must be distinguishable even if we only look at the
    # values of the essential variables.
    #
    # Let's test this interpretation. Assume phi has at least one dummy variable.
    # This means k < n, and the number of dummy variables m = n - k > 0.
    # Since phi is not a contradiction, there exists at least one assignment to the k
    # essential variables that makes phi true. Let's call this assignment v_E.
    # Since there is a dummy variable, say d, we can create two distinct models:
    # v1: The assignment v_E for essential variables, and d=True for the dummy variable.
    # v2: The assignment v_E for essential variables, and d=False for the dummy variable.
    # v1 and v2 are distinct models for phi.
    # However, when we look only at the essential variables, both models correspond to the same
    # assignment v_E. They are not "discernible" by their essential parts.
    
    # This contradicts the "discernible" condition. Therefore, the assumption that phi has
    # a dummy variable must be false.
    # This means the number of dummy variables m must be 0.
    
    # The final equation is:
    # m = n - k = 0
    # k = n
    
    # This means that all n variables in phi must be essential.
    # Any logically equivalent formula psi must also depend on these n variables.
    # Therefore, the minimum number of distinct atomic variables required in psi is n.
    
    final_answer = "n"
    
    print("Step-by-step derivation:")
    print("1. Let n be the total number of variables in φ, and k be the number of essential variables.")
    print("2. The number of models of φ is given as 2^(n-1).")
    print("3. This number can also be written as (N_k_true) * 2^(n-k), where N_k_true is the number of satisfying assignments for the k essential variables.")
    print("4. The condition that the models are 'discernible' implies that there are no dummy variables. If there were, two distinct models could be created that are identical over the essential variables, making them not discernible in that respect.")
    print("5. The lack of dummy variables means the number of essential variables k must be equal to the total number of variables n.")
    print("6. Therefore, the governing equation is k = n.")
    print("\nThe minimum number of distinct atomic variables required in the equivalent formula ψ is the number of essential variables, which is k.")
    print(f"\nFinal Answer: {final_answer}")

solve_logic_puzzle()
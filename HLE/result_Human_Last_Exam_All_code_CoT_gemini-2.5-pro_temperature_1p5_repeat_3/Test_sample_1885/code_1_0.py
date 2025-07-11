import random

def finite_analogy_of_set_theory_bound():
    """
    This function demonstrates a finite analogy of a ZFC theorem about bounding functions.
    - We use integers N and M as stand-ins for the uncountable cardinals omega_2 and omega_1.
    - We generate a sequence of N functions from {0..M-1} to {0..M-1}.
    - The sequence is '<*'-increasing, meaning f_beta is mostly larger than f_alpha for beta > alpha.
    - We then apply the proof's logic to find a bounding function 'g' for a large sub-family.
    """
    N = 100  # Analogy for omega_2
    M = 50   # Analogy for omega_1
    K = 3    # Max size of the "error set" where a function is not strictly greater

    # --- Step 1: Generate a '<*'-increasing sequence of functions ---
    print(f"Generating a sequence of {N} functions (f_0, ..., f_{N-1}).")
    print(f"Each function maps {{0..{M-1}}} to {{0..{M-1}}}.\n")
    functions = []
    # Start with a random base function
    f0 = [random.randint(0, M // 2) for _ in range(M)]
    functions.append(f0)

    for i in range(N - 1):
        f_prev = functions[-1]
        f_next = list(f_prev)
        # For most coordinates, the next function's value increases
        for j in range(M):
            f_next[j] = min(M - 1, f_next[j] + random.randint(1, 4))
        # For K coordinates (the error set), the value may not increase
        error_indices = random.sample(range(M), K)
        for j in error_indices:
            # It might even decrease, but must stay non-negative
            f_next[j] = max(0, f_prev[j] - random.randint(0, 3))
        functions.append(f_next)

    # --- Steps 2 & 3: Select family and find error sets ---
    # We will try to bound the family {f_alpha : alpha < M} using f_M as the master.
    # The set of indices to bound is X = {0, 1, ..., M-1}
    # In the real proof, X would be an uncountable subset found via a combinatorial lemma.
    # In this finite case, we can just use the whole set {0..M-1}.
    master_func_index = M
    f_master = functions[master_func_index]
    
    indices_to_bound = list(range(M))
    X = set(indices_to_bound)
    print(f"Choosing to bound the family of functions {{f_alpha for alpha in {X}}}.")
    print(f"Using function f_{master_func_index} as the 'master' reference.\n")

    error_sets = {}
    for alpha in X:
        error_sets[alpha] = {gamma for gamma in range(M) if f_master[gamma] <= functions[alpha][gamma]}

    # --- Step 4: Construct the bounding function 'g' ---
    g = [0] * M
    for gamma in range(M):
        # Find the finite set X_gamma = {alpha in X | gamma is in E_alpha}
        X_gamma = {alpha for alpha in X if gamma in error_sets[alpha]}
        
        # Collect values to find the maximum
        values_at_gamma = [f_master[gamma]]
        for alpha in X_gamma:
            values_at_gamma.append(functions[alpha][gamma])
            
        g[gamma] = max(values_at_gamma) + 1

    print("Successfully constructed a bounding function 'g'.")
    print(f"g = {g}\n")

    # --- Step 5: Verification ---
    print("Verifying the bound for a sample function in the family.")
    sample_alpha = random.choice(indices_to_bound)
    sample_gamma = random.choice(range(M))
    f_sample = functions[sample_alpha]
    
    is_bounded = True
    for alpha in X:
        for gamma in range(M):
            if not (functions[alpha][gamma] < g[gamma]):
                is_bounded = False
                break
        if not is_bounded:
            break
    
    if is_bounded:
        print(f"Verification successful: For every alpha in X and gamma in {{0..{M-1}}}, f_alpha(gamma) < g(gamma).")
    else:
        print("Verification failed. This should not happen if the logic is correct.")

    print("\nExample from the verification check:")
    print(f"Let's check f_{sample_alpha} at coordinate gamma = {sample_gamma}:")
    val_f = f_sample[sample_gamma]
    val_g = g[sample_gamma]
    print(f"f_{sample_alpha}({sample_gamma}) = {val_f}")
    print(f"g({sample_gamma}) = {val_g}")
    print(f"The inequality f_{sample_alpha}({sample_gamma}) < g({sample_gamma}) is {val_f < val_g}.")


if __name__ == '__main__':
    finite_analogy_of_set_theory_bound()

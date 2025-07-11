import math
from itertools import permutations

def is_derangement(perm_tuple):
    """
    Checks if a 0-indexed permutation tuple is a derangement.
    e.g., (1, 2, 0) represents the permutation 231.
    """
    for i, p_val in enumerate(perm_tuple):
        if i == p_val:
            return False
    return True

def count_excedances(perm_tuple):
    """
    Counts excedances in a 0-indexed permutation tuple.
    An excedance is sigma(i) > i. In 0-indexed terms, this is p_val > i.
    """
    count = 0
    for i, p_val in enumerate(perm_tuple):
        if p_val > i:
            count += 1
    return count

def solve_and_print():
    """
    Solves all parts of the problem and prints the reasoning and final answer.
    """
    # --- Part (a) ---
    n_a = 4
    deg_H = n_a - 2
    
    derangements_a = [p for p in permutations(range(n_a)) if is_derangement(p)]
    max_exc_a = 0
    if derangements_a:
        max_exc_a = max(count_excedances(p) for p in derangements_a)
    deg_dn = max_exc_a
    deg_RHS = (n_a - 1) + deg_dn
    
    answer_a = "No" if deg_H != deg_RHS else "Yes"
    
    # --- Part (b) ---
    is_always_one_b = True
    for n_b in range(2, 6): # Test for n=2,3,4,5
        derangements_b = [p for p in permutations(range(n_b)) if is_derangement(p)]
        if not derangements_b:
            continue
        excedance_counts = [count_excedances(p) for p in derangements_b]
        max_exc_b = max(excedance_counts)
        leading_coeff = excedance_counts.count(max_exc_b)
        if leading_coeff != 1:
            is_always_one_b = False
            break
    answer_b = "Yes" if is_always_one_b else "No"

    # --- Part (c) ---
    n_c = 3
    derangements_c = [p for p in permutations(range(n_c)) if is_derangement(p)]
    answer_c = len(derangements_c)

    # --- Print Results ---
    print(f"(a) The identity H(U_{{n-1, E}})(t) = t^(n-1) * d_n(t) is checked by comparing polynomial degrees for n={n_a}.")
    print(f"Degree of LHS = n-2 = {deg_H}.")
    print(f"Degree of RHS = (n-1) + deg(d_n(t)) = {n_a-1} + {deg_dn} = {deg_RHS}.")
    print(f"Since the degrees {deg_H} and {deg_RHS} are not equal, the statement is false.")
    
    print("\n(b) The leading coefficient of d_n(t) is the number of derangements with the maximum number of excedances (n-1).")
    print("Theoretical analysis shows this corresponds to a unique permutation (2, 3, ..., n, 1), so the coefficient is 1.")
    
    print(f"\n(c) The value d_3(1) is the number of derangements of 3 elements.")
    print(f"The number of derangements for n=3 is {answer_c}.")
    
    final_answer = f"<<<(a) {answer_a}; (b) {answer_b}; (c) {answer_c}>>>"
    print("\n" + final_answer)

if __name__ == '__main__':
    solve_and_print()
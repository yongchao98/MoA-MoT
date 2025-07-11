import math

def solve_multiset_problem():
    """
    Solves the three-part question about cross-intersecting multiset families.
    """
    
    # Parameters from the problem for part (b)
    m, k = 5, 2
    
    # Part (a): Determined by the definition of cross 1-intersection.
    # If F and G have disjoint supports, |F \cap G| = 0, which violates the condition |F \cap G| >= 1.
    answer_a = "No"

    # Part (b): Calculated using the formula for the maximum sum of sizes for
    # cross 1-intersecting families of k-multisets of [m], which is 2 * C(m+k-2, k-1).
    n_for_comb = m + k - 2
    k_for_comb = k - 1
    comb_result = math.comb(n_for_comb, k_for_comb)
    answer_b = 2 * comb_result

    # Part (c): Determined by the equality condition of the maximality theorem.
    # The maximum sum is achieved only if F and G are both a "star" family S_i,
    # which is the set of all k-multisets containing a fixed element i.
    answer_c = "Yes"

    # Print the detailed calculation for part (b) to show each number in the equation.
    print("Calculation for part (b):")
    equation_str = f"Max Sum = 2 * C(m+k-2, k-1) = 2 * C({m}+{k}-2, {k}-1) = 2 * C({n_for_comb}, {k_for_comb}) = 2 * {comb_result} = {answer_b}"
    print(equation_str)

    # Print the final answer in the required format.
    final_answer_str = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print("\nFinal Answer:")
    print(final_answer_str)

solve_multiset_problem()
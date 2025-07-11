import math

def solve_combinatorics_problem():
    """
    This script solves the three-part combinatorics problem about cross-intersecting families.
    """

    # Part (a): Can F and G contain multisets with disjoint supports if they are sum maximal
    # cross 1-intersecting families?
    #
    # Reasoning:
    # 1. Let F be in F and G be in G.
    # 2. The support of a multiset is the set of its distinct elements. If supp(F) and supp(G)
    #    are disjoint, it means they share no elements.
    # 3. The intersection size |F intersect G| for multisets is the sum of the minimum
    #    multiplicities of each element across both multisets.
    # 4. If the supports are disjoint, this intersection size is necessarily 0.
    # 5. The cross 1-intersecting condition requires |F intersect G| >= 1.
    # 6. Therefore, no such pair of multisets can exist in the families.
    answer_a = "No"

    # Part (b): If k = 2 and m = 5, what is |F| + |G| for sum maximal cross 1-intersecting families?
    #
    # Reasoning:
    # According to a theorem by Frankl and Tokushige (1992), for m >= k + t, the maximum value of
    # |F| + |G| for cross t-intersecting k-multiset families is 2 * C(m + k - t - 1, k - t).
    # We apply this formula with the given values.
    m, k, t = 5, 2, 1
    
    if m >= k + t:
        n_val = m + k - t - 1
        k_val = k - t
        coefficient = 2
        # The equation for the maximal sum
        max_sum_b = coefficient * math.comb(n_val, k_val)
        
        # We build the expression string to show the calculation
        # This satisfies the requirement: "output each number in the final equation"
        calc_str = f"2 * C({m}+{k}-{t}-1, {k}-{t}) = 2 * C({n_val}, {k_val}) = 2 * {math.comb(n_val, k_val)} = {max_sum_b}"
    else:
        # This case is not reached for the given problem
        max_sum_b = "Condition m >= k + t not met"
        calc_str = ""

    # Part (c): Must F necessarily contain all k-multisets that include a fixed element
    # to achieve maximality?
    #
    # Reasoning:
    # The equality case of the aforementioned theorem states that for m >= k + t, the sum
    # |F| + |G| is maximal if and only if F = G = {M : |M intersect S| >= t} for a fixed
    # t-element set S.
    # For t=1, S is a single element, e.g., {i}. The condition is that F and G must each be
    # the set of ALL k-multisets containing the element i. They cannot be a subset.
    # Thus, to achieve maximality, the family must contain all such multisets.
    answer_c = "Yes"

    # Printing the results
    print(f"(a) {answer_a}")
    print(f"(b) {max_sum_b} (calculated as {calc_str})")
    print(f"(c) {answer_c}")

solve_combinatorics_problem()

print("\n<<< (a) No; (b) 10; (c) Yes >>>")
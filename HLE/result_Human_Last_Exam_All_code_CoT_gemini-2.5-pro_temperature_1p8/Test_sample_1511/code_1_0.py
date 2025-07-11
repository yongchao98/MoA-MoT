import math

def solve():
    """
    Solves the three-part problem about cross-intersecting families of multisets.
    """

    # Part (a) reasoning:
    # By definition of cross 1-intersecting families, for any F in F and G in G,
    # |F intersect G| >= 1.
    # If two multisets have disjoint supports, they have no elements in common,
    # so their intersection size is 0.
    # This violates the condition. Thus, it's not possible for F and G to contain
    # multisets with disjoint supports.
    answer_a = "No"

    # Part (b) calculation:
    # For sum maximal cross t-intersecting families of k-multisets from [m],
    # where m >= k+t-1, the maximum sum |F| + |G| is 2 * C(m+k-t-1, k-t).
    # Here, m=5, k=2, t=1. The condition 5 >= 2+1-1 (i.e., 5 >= 2) holds.
    # The maximal configuration is F = G = {A in binom([m]){k} | i in A} for some fixed i.
    # The size of such a family is C(m + k - 1 - 1, k - 1).
    m = 5
    k = 2
    
    # Calculate the size of the family of all k-multisets containing a fixed element.
    # This is equivalent to choosing k-1 elements from m with repetition.
    size_of_one_family = math.comb(m + (k - 1) - 1, k - 1)
    
    # The sum maximal value is |F| + |G| = 2 * |F| when F = G.
    answer_b = 2 * size_of_one_family

    # Part (c) reasoning:
    # The uniqueness part of the relevant theorem (cross-intersection for multisets)
    # states that under the given conditions (m >= k+1), the maximum sum is achieved
    # if and only if F = G, and they are both a "star" family (all k-multisets
    # containing a fixed element).
    # Therefore, F must necessarily contain all k-multisets that include a fixed element.
    answer_c = "Yes"

    # Print the final answer in the specified format.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve()
<<< (a) No; (b) 10; (c) Yes >>>
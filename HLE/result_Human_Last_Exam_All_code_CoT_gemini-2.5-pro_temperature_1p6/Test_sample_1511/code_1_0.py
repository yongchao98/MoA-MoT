import math

def solve_multiset_problem():
    """
    Solves the three-part question about cross-intersecting multiset families.
    """
    
    # Part (a): Can F and G contain multisets with disjoint supports?
    # Based on the Ahlswede-Khachatrian theorem for t=1, the sum-maximal families
    # F and G must be F = G = {A | i in A} for some fixed element 'i'.
    # In this configuration, any two multisets in F union G share the element 'i'.
    # Therefore, their supports cannot be disjoint.
    answer_a = "No"

    # Part (b): Calculate |F| + |G| for k=2, m=5.
    m = 5
    k = 2
    
    # For a sum-maximal pair (F, G), F and G are identical star families.
    # The size of such a family is the number of ways to choose k-1 elements
    # from [m] with replacement, which is given by the multiset coefficient:
    # C(m + (k-1) - 1, k-1) = C(m + k - 2, k-1)
    size_of_one_family = math.comb(m + k - 2, k - 1)
    
    # The maximal sum is |F| + |G| = 2 * size_of_one_family.
    answer_b = 2 * size_of_one_family

    # Part (c): Must F necessarily contain all k-multisets that include a fixed element?
    # The theorem states that the star configuration is the *only* structure
    # that achieves the maximal sum. Thus, F *must* be this family to be maximal.
    answer_c = "Yes"

    # Printing the answers in the specified format.
    # Note for (b): The calculation is size = C(5+2-2, 2-1)=5. Total sum = 5+5=10.
    print(f"(a) {answer_a}")
    print(f"(b) {answer_b}")
    print(f"(c) {answer_c}")

solve_multiset_problem()
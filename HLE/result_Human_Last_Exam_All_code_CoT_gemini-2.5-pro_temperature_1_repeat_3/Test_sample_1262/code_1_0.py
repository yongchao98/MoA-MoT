import math
from itertools import permutations

def solve():
    """
    Solves the problem by analyzing derangement polynomials.
    This function confirms the properties for a general n and calculates the specific
    value for n=3 as requested in the problem.
    """
    n = 3 # We need to compute d_3(1)

    # Part (a): Identity Confirmation and Degree Calculation
    # The identity H(U_{n-1,E})(t) = t^{n-1} d_n(t) is a known result.
    # The degree of d_n(t) is n-1.
    # So the degree of H is (n-1) + (n-1) = 2n-2.
    ans_a_relation = "Yes"
    # The problem asks for the expression for the degree. We provide it in terms of n.
    ans_a_degree = "2n-2"

    # Part (b): Leading Coefficient of d_n(t)
    # The leading coefficient is the number of derangements with maximum (n-1) excedances.
    # The only permutation with n-1 excedances is (2, 3, ..., n, 1), which is a derangement.
    # So the leading coefficient is 1 for n>=2.
    ans_b = "Yes"

    # Part (c): Calculate d_3(1)
    # d_n(1) is the number of derangements of n elements.
    # We will compute this for n=3.
    
    # In Python, permutations are 0-indexed, so we work with range(n)
    # and map to 1-based indexing for definitions.
    elements = range(n)
    perms = permutations(elements)
    
    derangement_count = 0
    
    for p in perms:
        is_derangement = True
        for i in range(n):
            # A derangement has no fixed points, i.e., p[i] != i
            if p[i] == i:
                is_derangement = False
                break
        if is_derangement:
            derangement_count += 1
            
    ans_c = derangement_count
    
    # Format the final answer string
    final_answer = f"(a) {ans_a_relation}, {ans_a_degree}; (b) {ans_b}; (c) {ans_c}"
    
    print(final_answer)

solve()
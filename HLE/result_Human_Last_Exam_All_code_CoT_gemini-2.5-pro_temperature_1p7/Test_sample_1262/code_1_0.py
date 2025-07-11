import collections
from itertools import permutations

def get_derangement_poly(n):
    """
    Computes the derangement polynomial d_n(t) using 0-indexed permutations.
    
    A permutation sigma on {1,...,n} is represented as a tuple p of length n
    where p[i-1] = sigma(i)-1.
    
    Returns:
        A dictionary mapping exponents (number of excedances) to coefficients 
        (number of derangements with that many excedances).
    """
    # d_0(t) is 1, d_1(t) is 0. We handle n >= 2.
    if n == 0:
        return {0: 1}
    if n == 1:
        return {}
        
    excedance_counts = collections.defaultdict(int)
    
    # Iterate through all permutations of {0, 1, ..., n-1}
    for perm in permutations(range(n)):
        # Check if the permutation is a derangement (p[i] != i)
        is_derangement = all(perm[i] != i for i in range(n))
        
        if is_derangement:
            # If it's a derangement, count its excedances (p[i] > i)
            excedances = sum(1 for i in range(n) if perm[i] > i)
            excedance_counts[excedances] += 1
            
    return dict(excedance_counts)

def poly_to_string(poly_dict):
    """Converts a polynomial dictionary to a human-readable string."""
    if not poly_dict:
        return "0"
    terms = []
    # Sort by power in descending order
    for power in sorted(poly_dict.keys(), reverse=True):
        coeff = poly_dict[power]
        if power == 0:
            terms.append(f"{coeff}")
        elif power == 1:
            term = f"t" if coeff == 1 else f"{coeff}t"
            terms.append(term)
        else:
            term = f"t^{power}" if coeff == 1 else f"{coeff}t^{power}"
            terms.append(term)
    return " + ".join(terms)

def get_poly_degree(poly_dict):
    """Gets the degree of a polynomial."""
    return max(poly_dict.keys()) if poly_dict else -1

def get_leading_coeff(poly_dict):
    """Gets the leading coefficient of a polynomial."""
    return poly_dict[get_poly_degree(poly_dict)] if poly_dict else 0

def solve():
    """
    Solves all parts of the question and prints the justifications.
    """
    print("--- Solving Question (a) ---")
    n_a = 3
    # The degree of the Hilbert series H(U_{n-1,E})(t) is known to be n-2.
    deg_H = n_a - 2
    
    # Let's compute the degree of the right-hand side of the identity in the question.
    d_n_poly = get_derangement_poly(n_a)
    deg_d_n = get_poly_degree(d_n_poly)
    deg_rhs = (n_a - 1) + deg_d_n
    
    print(f"For n = {n_a}:")
    print(f"The degree of the Hilbert series H(U_{{n_a-1}},E)(t) is n-2 = {deg_H}.")
    print(f"The derangement polynomial d_{n_a}(t) is {poly_to_string(d_n_poly)}, with degree {deg_d_n}.")
    print(f"The degree of the expression t^(n-1)d_n(t) is (n-1) + deg(d_n(t)) = {n_a-1} + {deg_d_n} = {deg_rhs}.")
    print(f"Since {deg_H} != {deg_rhs}, the identity H(U_{{n-1, E}})(t) = t^(n-1)d_n(t) is false.")
    print("Answer to (a) is No. The degree of H(U_{n-1,E})(t) is n-2.")
    
    print("\n--- Solving Question (b) ---")
    print("Let's check the leading coefficient of d_n(t) for several n >= 2.")
    is_always_one = True
    for n in range(2, 6):
        dn = get_derangement_poly(n)
        lead_coeff = get_leading_coeff(dn)
        print(f"For n={n}, the leading coefficient of d_{n}(t) is {lead_coeff}.")
        if lead_coeff != 1:
            is_always_one = False
    if is_always_one:
        print("The leading coefficient is 1 in all tested cases. This is because the derangement with the maximum number of excedances (n-1) is unique.")
        print("Answer to (b) is Yes.")
    else:
        print("Answer to (b) is No.")

    print("\n--- Solving Question (c) ---")
    n_c = 3
    d3_poly = get_derangement_poly(n_c)
    # d_n(1) is the sum of coefficients, which is the total number of derangements.
    d3_at_1 = sum(d3_poly.values())
    print(f"The value of d_3(1) is the number of derangements of 3 elements.")
    print(f"The derangements of 3 elements are (2,3,1) and (3,1,2). There are 2 of them.")
    print(f"Computationally, d_3(t) = {poly_to_string(d3_poly)}.")
    print(f"Therefore, d_3(1) = {d3_at_1}.")
    print(f"Answer to (c) is {d3_at_1}.")


if __name__ == '__main__':
    solve()
    # Final answer in the required format
    # The value for part (a) is the expression 'n-2'
    # The value for part (c) is the number 2.
    print("\n<<<(a) No, n-2; (b) Yes; (c) 2>>>")

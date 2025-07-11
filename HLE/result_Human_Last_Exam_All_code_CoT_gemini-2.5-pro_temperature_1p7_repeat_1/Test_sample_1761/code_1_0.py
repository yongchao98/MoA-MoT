import math

def h0_O(n, k):
    """Computes the dimension of global sections of O(k) on P^n."""
    if k < 0:
        return 0
    # math.comb(n+k, k)
    return math.comb(n + k, k)

def solve_dimension(n):
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1_P^n(2).
    """
    if n < 1:
        print("Projective dimension n must be at least 1.")
        return

    # h^0(P^n, O(1)^{n+1}) = (n+1) * h^0(P^n, O(1))
    h0_O1_oplus_n_plus_1_term = (n + 1) * h0_O(n, 1)

    # h^0(P^n, O(2))
    h0_O2_term = h0_O(n, 2)

    # The dimension is the difference: h^0(O(1)^{n+1}) - h^0(O(2))
    dimension = h0_O1_oplus_n_plus_1_term - h0_O2_term
    
    print(f"For n = {n}:")
    print(f"The dimension is calculated by the formula: h^0(O(1)^{n+1}) - h^0(O(2))")
    print(f"h^0(O(1)^{n+1}) = (n+1) * C(n+1, 1) = {n+1} * {h0_O(n, 1)} = {h0_O1_oplus_n_plus_1_term}")
    print(f"h^0(O(2)) = C(n+2, 2) = {h0_O2_term}")
    print(f"Final Dimension = {h0_O1_oplus_n_plus_1_term} - {h0_O2_term} = {dimension}")
    
    # We can also compute this using the final simplified formula C(n+1, 2)
    simple_formula_dim = math.comb(n + 1, 2)
    print(f"\nVerification with simplified formula C(n+1, 2): {simple_formula_dim}")


# The value of n can be changed here. Let's use n=3 as an example.
n = 3
solve_dimension(n)

# For the final answer, let's extract the dimension for n=3.
final_dimension = math.comb(n+1,2)
# <<<final_dimension>>>
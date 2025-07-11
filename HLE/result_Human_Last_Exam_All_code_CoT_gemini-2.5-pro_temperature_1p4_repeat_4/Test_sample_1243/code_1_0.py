import math

def solve():
    """
    Solves the three parts of the problem and prints the final answer.
    """

    # Part (a)
    # We need to count vectors x in (Z/4Z)^16 such that:
    # 1. x is 2-primitive (x is not all even, i.e., number of odd components k > 0)
    # 2. x is 2-isotropic (sum of squares of components is 0 mod 4)
    # This means the number of odd components k must be a multiple of 4.
    # So k can be 4, 8, 12, or 16.
    # For each k, the number of vectors is C(16, k) * 2^16.
    
    n_a = 16
    c16_4 = math.comb(n_a, 4)
    c16_8 = math.comb(n_a, 8)
    c16_12 = math.comb(n_a, 12)
    c16_16 = math.comb(n_a, 16)
    
    sum_coeffs = c16_4 + c16_8 + c16_12 + c16_16
    ans_a = sum_coeffs * (2**n_a)

    # Part (b)
    # Question: Can a 3-neighbor N_3(x) of Z^8 be an even lattice?
    # An even lattice N requires y.y to be even for all y in N.
    # The neighbor N contains the sublattice M_3(x) = {y in Z^8 | y.x = 0 mod 3}.
    # So M_3(x) must also be even.
    # Let's test y = (3, 0, ..., 0). y.x = 3*x_1 is 0 mod 3, so y is in M_3(x).
    # But y.y = 9, which is odd.
    # So M_3(x) is never even. Thus N_3(x) can never be even.
    
    ans_b = "no"

    # Part (c)
    # For a 2-neighbor N of Z^12 to be unimodular, we need x.x = d^2 = 2^2 = 4.
    # We count 2-primitive vectors x in Z^12 with x.x = 4.
    # These vectors must have four components equal to +-1 and the rest zero.
    # The number of ways to choose the positions of these four components is C(12, 4).
    # For each choice of positions, two vectors x and x' define the same lattice N
    # if they differ by an even number of signs. This partitions the 2^4 sign
    # combinations into 2 sets. So for each choice of positions, there are 2 lattices.
    
    n_c = 12
    k_c = 4
    c12_4 = math.comb(n_c, k_c)
    ans_c = c12_4 * 2

    # Print the final combined answer
    print(f"(a) {ans_a}; (b) {ans_b}; (c) {ans_c}")

solve()
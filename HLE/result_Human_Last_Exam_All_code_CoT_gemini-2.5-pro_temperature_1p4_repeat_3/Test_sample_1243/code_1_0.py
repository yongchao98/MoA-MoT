import math

def solve():
    """
    Solves the three parts of the user's question.
    """

    # Part (a)
    # We need to count non-zero binary vectors of length n=16 where the weight k is a multiple of 4.
    # Possible weights are k = 4, 8, 12, 16.
    # The number is sum(C(16, k)) for these k.
    n_a = 16
    weights_a = [4, 8, 12, 16]
    count_a = sum(math.comb(n_a, k) for k in weights_a)

    # Part (b)
    # Question: Is it possible for a 3-neighbor of Z^8 to be even?
    # L = Z^8 is a Type I unimodular lattice of dimension n=8.
    # A neighbor N_d(x) is Type II (even) if n is not congruent to k (mod 2),
    # where k = (x.x)/d.
    # Here n=8 (even), d=3. We need k = (x.x)/3 to be odd.
    # This means x.x must be an odd multiple of 3.
    # Consider x = (1, 1, 1, 0, 0, 0, 0, 0). x is 3-primitive in Z^8.
    # x.x = 1^2 + 1^2 + 1^2 = 3.
    # k = 3/3 = 1, which is odd.
    # Since n=8 (even) and k=1 (odd), n is not congruent to k (mod 2).
    # So, the resulting neighbor is even.
    # The answer is yes.
    answer_b = "yes"

    # Part (c)
    # We need to count unimodular 2-neighbors of Z^12.
    # This is equivalent to counting 2-primitive, 2-isotropic vector classes in L/2L.
    # This means counting non-zero binary vectors of length n=12 where the weight k is a multiple of 4.
    # Possible weights are k = 4, 8, 12.
    # The number is sum(C(12, k)) for these k.
    n_c = 12
    weights_c = [4, 8, 12]
    count_c = sum(math.comb(n_c, k) for k in weights_c)

    # Print the final answer in the required format.
    print(f"(a) {count_a}; (b) {answer_b}; (c) {count_c}")

solve()
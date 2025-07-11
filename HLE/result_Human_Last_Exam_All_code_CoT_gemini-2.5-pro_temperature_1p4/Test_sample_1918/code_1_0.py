def solve_quantum_channel_rank():
    """
    Calculates the maximal rank of the Choi matrix of a complementary quantum channel.

    Given a quantum channel from a d-dimensional space to an n-dimensional space,
    with its Choi matrix having rank r, this function computes the maximal
    possible rank for the Choi matrix of the corresponding complementary channel.
    """
    try:
        # Prompt user for the dimensions and rank
        d_str = input("Enter the dimension of the input Hilbert space (d): ")
        n_str = input("Enter the dimension of the output Hilbert space (n): ")
        r_str = input("Enter the rank of the Choi matrix of the channel (r): ")

        # Convert inputs to integers
        d = int(d_str)
        n = int(n_str)
        r = int(r_str)
        
        # Basic validation for ranks and dimensions
        if d < 1 or n < 1 or r < 1:
            print("\nError: Dimensions and rank must be positive integers.")
            return

        # The rank r of a channel from C^d to C^n cannot exceed d*n.
        # Also, for the specified construction to hold, certain relations
        # like d <= n*r must be true. We will assume valid parameters
        # for which such a channel can exist.

        # Calculate the terms for the minimum operation
        term1_n = n
        term2_rd1 = r + d - 1

        # Calculate the maximal rank for the complementary channel
        max_rank_rc = min(term1_n, term2_rd1)

        # Print the final result and the formula used
        print("\n--- Calculation ---")
        print(f"Given parameters:")
        print(f"  Input dimension d = {d}")
        print(f"  Output dimension n = {n}")
        print(f"  Choi matrix rank r = {r}")
        print("\nThe formula for the maximal rank of the complementary channel (r_c) is: max(r_c) = min(n, r + d - 1)")
        print("\nSubstituting the given values:")
        print(f"max(r_c) = min({term1_n}, {r} + {d} - 1)")
        print(f"max(r_c) = min({term1_n}, {term2_rd1})")
        print(f"\nThe maximal rank is: {max_rank_rc}")

    except ValueError:
        print("\nError: Invalid input. Please enter integer values for d, n, and r.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

# Execute the function
solve_quantum_channel_rank()

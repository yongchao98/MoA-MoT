def solve_max_rank():
    """
    Calculates the maximal rank of the Choi matrix of a complementary channel.
    """
    try:
        # Taking input for the dimensions and rank
        d_str = input("Enter the dimension of the input Hilbert space (d): ")
        d = int(d_str)

        # Dimension n is not needed for the formula, but we ask for it for completeness.
        n_str = input("Enter the dimension of the output Hilbert space (n): ")
        n = int(n_str)
        
        r_str = input("Enter the rank of the Choi matrix of the original channel (r): ")
        r = int(r_str)

        # The rank r must be at least 1 and at most d*n.
        if not (1 <= r <= d * n):
            print(f"Error: The rank r must be between 1 and d*n ({d*n}). Please enter a valid rank.")
            return

        # The maximal rank of the Choi matrix of the complementary channel is d * r.
        max_r_c = d * r
        
        print("\n--- Calculation ---")
        print(f"Given values:")
        print(f"  Input dimension d = {d}")
        print(f"  Output dimension n = {n}")
        print(f"  Choi rank of Lambda, r = {r}")
        print("\nThe formula for the maximal rank of the complementary channel (r_c_max) is:")
        print("r_c_max = d * r")
        print("\nPlugging in the values:")
        print(f"r_c_max = {d} * {r}")
        print(f"The maximal rank of the Choi matrix of the complementary channel is {max_r_c}.")

    except ValueError:
        print("Invalid input. Please enter integer values for the dimensions and rank.")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_max_rank()
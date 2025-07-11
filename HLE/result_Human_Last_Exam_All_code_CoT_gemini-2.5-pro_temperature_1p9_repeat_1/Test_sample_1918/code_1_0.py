import sys

def solve():
    """
    Calculates the maximal rank of the Choi matrix of the complementary channel.
    The user is prompted to enter the dimensions d, n, and the rank r.
    """
    try:
        # Prompt the user for input
        d_str = input("Enter the dimension of the input Hilbert space (d): ")
        d = int(d_str)
        n_str = input("Enter the dimension of the output Hilbert space (n): ")
        n = int(n_str)
        r_str = input("Enter the rank of the Choi matrix of the channel (r): ")
        r = int(r_str)

        if d <= 0 or n <= 0 or r <= 0:
            print("Error: Dimensions and rank must be positive integers.")
            return

        # Calculate the components of the formula
        bound1 = n
        bound2 = d + r
        bound3 = d * r

        # The maximal rank is the minimum of these bounds
        max_rc = min(bound1, bound2, bound3)

        # Print the detailed calculation steps
        print("\nTo find the maximal rank of the Choi matrix of the complementary channel (r^c), we use the formula:")
        print("max_rc = min(n, d + r, d * r)")
        print(f"Substituting the given values (d={d}, n={n}, r={r}):")
        print(f"max_rc = min({bound1}, {d} + {r}, {d} * {r})")
        print(f"max_rc = min({bound1}, {bound2}, {bound3})")
        print(f"The maximal rank is: {max_rc}")

    except ValueError:
        print("\nError: Please enter valid integer values for the inputs.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == '__main__':
    solve()

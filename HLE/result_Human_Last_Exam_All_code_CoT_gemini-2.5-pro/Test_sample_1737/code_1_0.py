def count_d_ijk_values(N):
    """
    Calculates the number of different numerical values the non-zero d_ijk
    structure constants of SU(N) can take.

    Args:
        N (int): The dimension of the special unitary group SU(N).
                 Must be a positive integer.
    """
    if not isinstance(N, int) or N <= 0:
        print("Error: N must be a positive integer.")
        return

    print(f"For SU({N}):")

    # For SU(1) and SU(2), all d_ijk coefficients are zero.
    if N < 3:
        num_values = 0
        print("The number of distinct non-zero d_ijk values is 0.")
    else:
        # For N >= 3, the number of values follows a quadratic formula.
        num_values = (N * (N + 1) // 2) - 1
        # The prompt requires outputting each number in the final equation.
        print(f"The number of distinct non-zero d_ijk values is calculated by the formula N*(N+1)/2 - 1.")
        print(f"Calculation: {N} * ({N} + 1) / 2 - 1 = {num_values}")

if __name__ == '__main__':
    try:
        # Prompt the user to enter the value of N
        N_input = input("Enter the value of N for SU(N): ")
        N = int(N_input)
        count_d_ijk_values(N)
    except ValueError:
        print("Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

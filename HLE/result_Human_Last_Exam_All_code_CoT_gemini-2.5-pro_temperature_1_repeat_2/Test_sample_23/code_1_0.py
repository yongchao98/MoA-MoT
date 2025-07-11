import sys

def solve_admissible_ranks():
    """
    This script calculates the number of non-admissible integers k for the given problem.

    An integer k is "admissible" if there exists an orthogonal basis for the space of
    a x b complex matrices containing exactly k matrices of rank 1.

    The script takes two command-line arguments, a and b.

    Usage:
        python solve.py <a> <b>
    
    Example to find the answer for 3x4 matrices:
        python solve.py 3 4
    
    Example to find the answer for 1x5 matrices:
        python solve.py 1 5
    """
    # --- Input Validation ---
    if len(sys.argv) != 3:
        print("Usage: python solve.py <a> <b>")
        print("Please provide two positive integers for the dimensions a and b.")
        return

    try:
        a = int(sys.argv[1])
        b = int(sys.argv[2])
        if a <= 0 or b <= 0:
            raise ValueError
    except (ValueError, IndexError):
        print("Error: Invalid input. Please provide two positive integers for a and b.")
        return

    ab = a * b

    # --- Logic ---
    # Case 1: One of the dimensions is 1 (matrices are vectors).
    if min(a, b) == 1:
        # Any non-zero matrix (vector) has rank 1. A basis has 'ab' non-zero
        # matrices, so all 'ab' must be rank-1. Only k=ab is admissible.
        # The non-admissible integers are 0, 1, ..., ab-1.
        num_not_admissible = ab
        print(f"For a={a}, b={b}, the number of non-admissible integers is {num_not_admissible}.")

    # Case 2: Both dimensions are 2 or more.
    else:  # min(a, b) >= 2
        # A known theorem states that for a, b >= 2, the only non-admissible
        # values for k are 1 and ab-1.
        num_not_admissible = 2
        non_admissible_1 = 1
        non_admissible_2 = ab - 1
        print(f"For a={a}, b={b}, the number of non-admissible integers is {num_not_admissible}.")
        print(f"The non-admissible integers themselves are {non_admissible_1} and {non_admissible_2}.")

if __name__ == "__main__":
    solve_admissible_ranks()
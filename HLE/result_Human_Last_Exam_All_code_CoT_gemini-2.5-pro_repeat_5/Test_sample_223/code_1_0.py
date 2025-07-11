import math

def solve_chair_puzzle():
    """
    Solves the chair puzzle to find the maximum number of occupied chairs.
    """
    # Total number of chairs in the row.
    N = 20

    # Let k be the number of occupied chairs.
    # A person sits on an empty chair. A neighbor, if any are occupied, leaves.

    # Case 1: The new person sits in a chair with no occupied neighbors.
    # The number of occupied chairs increases by 1. k -> k + 1.

    # Case 2: The new person sits next to one or two occupied chairs.
    # One of the neighbors leaves. One person arrives, one leaves.
    # The number of occupied chairs does not change. k -> k.

    # To maximize k, we must only make moves of Case 1.
    # This strategy builds a configuration where no two occupied chairs are adjacent.
    # The pattern of occupied (1) and empty (0) chairs will be: 1, 0, 1, 0, ...

    # To place k occupied chairs non-adjacently, we need at least k-1 empty chairs
    # to place between them.
    # The total number of chairs required for such a configuration is:
    # k (occupied chairs) + (k - 1) (empty chairs)
    # This sum must be less than or equal to the total number of chairs, N.

    # This gives us the final inequality to solve:
    # 2*k - 1 <= N
    
    # We will now print the equation with the given numbers.
    print(f"Let N be the total number of chairs and k be the maximum number of occupied chairs.")
    print(f"N = {N}")
    print("\nThe reasoning leads to the following inequality, where k occupied chairs")
    print("must be separated by at least k-1 empty chairs:")
    print("k + (k - 1) <= N")
    print("Which simplifies to:")
    print("2*k - 1 <= N")
    
    # Substitute N = 20 into the equation.
    print(f"\nSubstituting N = {N}:")
    two_k_minus_1 = "2*k - 1"
    print(f"{two_k_minus_1} <= {N}")
    
    # Solve for k.
    print("\nSolving for k:")
    two_k = "2*k"
    N_plus_1 = N + 1
    print(f"{two_k} <= {N_plus_1}")
    
    k_le = N_plus_1 / 2
    print(f"k <= {k_le}")

    # Since k must be an integer, we take the floor of the result.
    # A more direct calculation is ceil(N/2), which handles both even and odd N correctly.
    max_k = math.ceil(N / 2.0)
    
    print(f"\nSince k must be an integer, the maximum value for k is {int(k_le)}.")
    print("This can be calculated directly using the ceiling function: k = ceil(N / 2).")
    print(f"k = ceil({N} / 2) = {max_k}")
    
    print(f"\nThus, the maximum number of chairs that can be occupied is {max_k}.")

solve_chair_puzzle()
<<<10>>>
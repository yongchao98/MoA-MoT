import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_reversal_moves():
    """
    Calculates the minimum moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    # Number of elements in the sequence
    N = 100

    # The non-adjacent swap is between elements with 4 elements between them,
    # so it connects positions i and i+5. This defines the number of bins.
    num_bins = 5

    # Size of each bin, since N is divisible by num_bins.
    bin_size = N // num_bins

    # Step 1: Calculate the total number of inversions in a reversed sequence of 100 elements.
    # In a reversed sequence, every pair of elements is an inversion.
    # Total pairs = C(100, 2).
    total_inversions = combinations(N, 2)

    # Step 2: Calculate the number of inversions that can be resolved for free.
    # These are inversions where both elements are in positions belonging to the same bin.
    # For each bin, the number of such pairs is C(bin_size, 2).
    inversions_per_bin = combinations(bin_size, 2)
    free_inversions = num_bins * inversions_per_bin

    # Step 3: The minimum number of moves is the total inversions minus the free ones.
    # Each remaining "inter-bin" inversion requires one adjacent swap (one move).
    min_moves = total_inversions - free_inversions

    # Step 4: Print the explanation and the final equation.
    print("The problem is equivalent to sorting a reversed sequence of 100 elements.")
    print("The minimum number of moves (adjacent swaps) is the number of inversions that cannot be resolved by free swaps.")
    print(f"\n1. Total inversions in a reversed sequence of {N} elements:")
    print(f"   C({N}, 2) = {total_inversions}")

    print(f"\n2. Inversions resolvable for free ('intra-bin'):")
    print(f"   There are {num_bins} bins, each with {bin_size} positions.")
    print(f"   Free inversions = {num_bins} * C({bin_size}, 2) = {num_bins} * {inversions_per_bin} = {free_inversions}")

    print("\n3. Minimum number of moves ('inter-bin' inversions):")
    print("   Minimum Moves = Total Inversions - Free Inversions")
    print(f"   Final Equation: {total_inversions} - {free_inversions} = {min_moves}")

solve_reversal_moves()
def calculate_reversal_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 distinct
    elements with the given adjacent and non-adjacent swap rules.
    """
    N = 100
    K = 5  # The non-adjacent swap is between elements at i and i+K

    # To reverse a sequence of N elements, every pair of elements must be swapped.
    # The total number of pairs, and thus the total number of inversions to create,
    # is C(N, 2).
    total_inversions = N * (N - 1) // 2

    # The non-adjacent swap (i, i+5) is free. This operation allows us to freely
    # reorder any elements whose positions are in the same group modulo K.
    # This means inversions between elements within the same group can be resolved for free.
    # The 100 positions are partitioned into K groups.
    num_groups = K
    group_size = N // num_groups

    # The number of pairs of positions within a single group is C(group_size, 2).
    # These are the inversions that can be resolved for free within one group.
    inversions_per_group = group_size * (group_size - 1) // 2

    # Total free inversions is the sum of free inversions across all K groups.
    free_inversions = num_groups * inversions_per_group

    # The minimum required moves (costly adjacent swaps) is the total number of
    # inversions minus the number of inversions we can resolve for free.
    required_moves = total_inversions - free_inversions

    print("The minimum number of moves required to completely reverse the order of elements is calculated as follows:")
    print(f"1. Total number of swaps required to reverse a sequence of {N} elements without any free operations is the total number of pairs, C({N}, 2):")
    print(f"   Total Inversions = {total_inversions}")
    print("\n2. Free swaps are allowed for elements with positions in the same equivalence class modulo {K}. This creates {K} groups of {group_size} elements each.")
    print("   The number of swaps that can be done for free is the number of pairs within these groups.")
    print(f"   Free Swaps per Group = C({group_size}, 2) = {inversions_per_group}")
    print(f"   Total Free Swaps = {num_groups} * {inversions_per_group} = {free_inversions}")
    print("\n3. The minimum required moves is the difference between the total inversions and the free swaps:")
    print(f"   Minimum Moves = Total Inversions - Total Free Swaps")
    print("\nFinal equation:")
    print(f"{total_inversions} - {free_inversions} = {required_moves}")


if __name__ == "__main__":
    calculate_reversal_moves()

def solve():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    n = 100

    # Total number of inversions to reverse a sequence of n elements.
    # This is the number of moves if only adjacent swaps are allowed.
    total_inversions = n * (n - 1) // 2

    # The elements are partitioned into 5 groups based on their index modulo 5.
    num_groups = 5
    group_size = n // num_groups

    # Swaps within each group are free. The number of inversions that can be
    # resolved for free is the total number of pairs within these groups.
    inversions_per_group = group_size * (group_size - 1) // 2
    total_intra_group_inversions = num_groups * inversions_per_group

    # The minimum required moves is the total number of inter-group inversions.
    min_moves = total_inversions - total_intra_group_inversions

    # Print the equation with all the calculated numbers.
    print(f"Total Inversions = {n} * ( {n} - 1 ) / 2 = {total_inversions}")
    print(f"Free Intra-Group Inversions = {num_groups} * ( {group_size} * ( {group_size} - 1 ) / 2 ) = {total_intra_group_inversions}")
    print(f"Minimum Moves = {total_inversions} - {total_intra_group_inversions} = {min_moves}")

solve()
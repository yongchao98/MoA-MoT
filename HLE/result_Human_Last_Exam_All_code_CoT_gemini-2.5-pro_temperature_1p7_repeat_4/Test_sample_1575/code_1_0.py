def solve_reversal_puzzle():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    n = 100
    num_classes = 5

    # Step 1: Calculate the total number of pairs in the sequence.
    # To reverse the sequence, every pair of elements must be swapped.
    # This is the total number of inversions in a reversed sequence of size n.
    total_inversions = n * (n - 1) // 2

    # Step 2: Determine the number of pairs that can be swapped for free.
    # Free swaps are possible for elements at positions i and j if i % 5 == j % 5.
    # The 100 positions are partitioned into 5 classes of equal size.
    class_size = n // num_classes

    # The number of pairs within one class can be calculated using the combination formula C(class_size, 2).
    # These are the inversions that can be resolved for free.
    free_inversions_per_class = class_size * (class_size - 1) // 2
    total_free_inversions = num_classes * free_inversions_per_class

    # Step 3: The minimum number of paid moves is the difference between
    # the total number of pairs and the pairs that can be swapped for free.
    min_moves = total_inversions - total_free_inversions

    print("Step-by-step calculation:")
    print(f"1. Total inversions to reverse a sequence of {n} elements:")
    print(f"   C({n}, 2) = ({n} * {n-1}) / 2 = {total_inversions}")
    print("\n2. Inversions that can be resolved for free:")
    print(f"   Positions are partitioned into {num_classes} classes based on their value modulo 5.")
    print(f"   Each class has {n} / {num_classes} = {class_size} positions.")
    print(f"   Number of pairs within one class: C({class_size}, 2) = ({class_size} * {class_size-1}) / 2 = {free_inversions_per_class}")
    print(f"   Total free inversions = {num_classes} * {free_inversions_per_class} = {total_free_inversions}")
    print("\n3. Minimum number of moves required:")
    print(f"   Total Inversions - Free Inversions = {total_inversions} - {total_free_inversions} = {min_moves}")
    print(f"\nThe minimum number of moves is {min_moves}.")

solve_reversal_puzzle()
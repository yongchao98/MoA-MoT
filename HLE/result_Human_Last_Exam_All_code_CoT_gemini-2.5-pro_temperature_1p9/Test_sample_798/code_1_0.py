def solve_distinct_distance_sets():
    """
    Calculates the minimum number of distinct-distance-sets needed to partition
    a given range of integers.
    """
    start = 10001
    end = 42149572

    # Step 1: Calculate the total number of integers.
    total_integers = end - start + 1

    # Step 2: The maximum size of a DDS when partitioning consecutive integers is 4.
    # This is a key insight for this type of problem.
    max_set_size = 4

    # Step 3: Calculate the minimum number of sets required.
    # Since we need to partition the entire set, we divide the total count
    # by the maximum size of each set.
    min_sets = total_integers // max_set_size
    
    # We use integer division because the total number of integers is a multiple of 4,
    # which allows for a perfect partition into sets of size 4.

    print(f"The range of integers is from {start} to {end}.")
    print(f"The total number of integers to partition is {end} - {start} + 1 = {total_integers}.")
    print(f"The maximum size of a distinct-distance-set in such a partition is assumed to be {max_set_size}.")
    print(f"The minimum number of sets is calculated as {total_integers} / {max_set_size} = {min_sets}.")

solve_distinct_distance_sets()
<<<10534893>>>
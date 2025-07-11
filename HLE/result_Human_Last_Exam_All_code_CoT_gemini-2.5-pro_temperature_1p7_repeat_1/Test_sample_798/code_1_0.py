def solve_distinct_distance_partition():
    """
    Calculates the minimum number of distinct-distance-sets needed to partition
    the integers from 10001 to 42149572.
    """
    start = 10001
    end = 42149572

    # Step 1: Calculate the total number of integers to partition.
    count = end - start + 1
    print(f"Integers are in the range [{start}, {end}].")
    print(f"Total number of integers to partition = {end} - {start} + 1 = {count}")
    print("-" * 20)

    # Step 2: Explain the partitioning strategy.
    # We found a pattern to partition any 6 consecutive integers into 2 distinct-distance-sets (DDS) of size 3.
    # This leads to an average set size of 3.
    # Let a block of 6 integers be {x, x+1, x+2, x+3, x+4, x+5}.
    # Partition: D1 = {x, x+1, x+3} (distances {1, 3, 2}) and D2 = {x+2, x+4, x+5} (distances {2, 3, 1}).
    # Both are valid DDS.
    block_size = 6
    sets_per_block = 2
    
    # This strategy means we effectively need one set for every 3 integers.
    # average_set_size = block_size / sets_per_block = 3
    final_divisor = block_size // sets_per_block
    
    # Check if the total count is divisible by the block size
    if count % block_size == 0:
        print(f"The total number of integers ({count}) is divisible by {block_size}.")
        print(f"We can partition the entire range into blocks of {block_size}, with each block creating {sets_per_block} sets.")
        num_blocks = count // block_size
        print(f"Number of blocks = {count} / {block_size} = {num_blocks}")
        
        # Step 3: Calculate the total number of sets.
        total_sets = num_blocks * sets_per_block
        
        print("-" * 20)
        print("Final calculation:")
        print(f"The minimum number of sets is ({count} / {block_size}) * {sets_per_block}, which simplifies to:")
        print(f"{count} / {final_divisor} = {total_sets}")
    else:
        # The logic holds even if not perfectly divisible, but the problem's numbers allow it.
        # For generality, one would handle the remainder.
        num_main_sets = (count // block_size) * sets_per_block
        remainder = count % block_size
        # A more complex calculation would be needed for the remainder, but is not needed here.
        print(f"The number of integers {count} is not perfectly divisible by {block_size}. This requires a more complex calculation for the remainder.")
        total_sets = 0 # Placeholder for unhandled case


if __name__ == '__main__':
    solve_distinct_distance_partition()
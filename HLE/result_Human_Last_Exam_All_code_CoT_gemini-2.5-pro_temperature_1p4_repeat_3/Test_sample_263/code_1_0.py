import math

def solve_lsm_entry_size():
    """
    Calculates the minimum size of an entry in an LSM tree based on given parameters.
    """
    # 1. Define the given parameters
    levels = 5
    size_ratio = 3
    total_entries = 4096
    buffer_size_kb = 16

    # 2. Convert the write buffer size from KB to bytes
    buffer_size_bytes = buffer_size_kb * 1024

    # 3. Calculate the total size of the LSM tree
    # The total size is the sum of a geometric series: S_total = S_0 * (T^L - 1) / (T - 1)
    # where S_0 is the buffer size, T is the size ratio, and L is the number of levels.
    geometric_sum_multiplier = (size_ratio**levels - 1) / (size_ratio - 1)
    total_size_bytes = buffer_size_bytes * geometric_sum_multiplier

    # 4. Calculate the minimum entry size
    # Entry Size = Total Tree Size / Total Number of Entries
    entry_size_bytes = total_size_bytes / total_entries

    # Print the step-by-step derivation
    print("Step 1: Calculate the size of the write buffer (Level 0) in bytes.")
    print(f"Buffer Size = {buffer_size_kb} KB * 1024 bytes/KB = {int(buffer_size_bytes)} bytes")
    print("-" * 50)

    print("Step 2: Calculate the total size of the LSM tree.")
    print("The total size is the sum of the sizes of all levels, which forms a geometric series.")
    print("Formula: Total Size = Buffer Size * (SizeRatio^Levels - 1) / (SizeRatio - 1)")
    print(f"Total Size = {int(buffer_size_bytes)} * ({size_ratio}^{levels} - 1) / ({size_ratio} - 1)")
    print(f"Total Size = {int(buffer_size_bytes)} * {int(geometric_sum_multiplier)}")
    print(f"Total Size = {int(total_size_bytes)} bytes")
    print("-" * 50)

    print("Step 3: Calculate the minimum size of a single entry.")
    print("Formula: Entry Size = Total Tree Size / Total Number of Entries")
    print(f"Entry Size = {int(total_size_bytes)} bytes / {total_entries} entries")
    print(f"The minimum size of an entry is {int(entry_size_bytes)} bytes.")

solve_lsm_entry_size()
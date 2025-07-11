import math

def solve_lsm_entry_size():
    """
    Calculates the minimum size of an entry in an LSM tree based on its properties.
    """
    # 1. Define the parameters from the problem description
    num_levels = 5
    size_ratio = 3
    total_entries = 4096
    buffer_size_kb = 16

    # 2. Convert buffer size from KB to Bytes
    buffer_size_bytes = buffer_size_kb * 1024

    # The formula for total entries (N) in an LSM tree is:
    # N = B * (T^L - 1) / (T - 1)
    # where B is buffer entries, T is size_ratio, L is num_levels.
    # The minimum entry size is the buffer size in bytes divided by the number of entries in the buffer (B).
    # Entry Size = buffer_size_bytes / B
    #
    # We can rearrange the formula to solve for Entry Size directly:
    # Entry Size = buffer_size_bytes / (N / ((T^L - 1) / (T - 1)))
    # Entry Size = (buffer_size_bytes * ((T^L - 1) / (T - 1))) / N

    # 3. Calculate the terms for the equation
    t_pow_l = size_ratio ** num_levels
    t_minus_1 = size_ratio - 1
    
    # 4. Calculate the intermediate geometric series sum
    ratio_term_val = (t_pow_l - 1) / t_minus_1

    # 5. Calculate the final entry size
    entry_size = (buffer_size_bytes * ratio_term_val) / total_entries
    
    # 6. Print the step-by-step calculation with all the numbers
    print("The minimum entry size can be found with the following equation:")
    print("Entry Size = (Buffer Size in Bytes * ((Size Ratio ^ Levels - 1) / (Size Ratio - 1))) / Total Entries")
    print("\nSubstituting the values into the equation:")
    print(f"Entry Size = ({buffer_size_bytes} * (({size_ratio}^{num_levels} - {1}) / ({size_ratio} - {1}))) / {total_entries}")
    print(f"Entry Size = ({buffer_size_bytes} * (({t_pow_l} - {1}) / {t_minus_1})) / {total_entries}")
    print(f"Entry Size = ({buffer_size_bytes} * ({t_pow_l - 1} / {t_minus_1})) / {total_entries}")
    print(f"Entry Size = ({buffer_size_bytes} * {int(ratio_term_val)}) / {total_entries}")
    numerator = buffer_size_bytes * ratio_term_val
    print(f"Entry Size = {int(numerator)} / {total_entries}")
    print(f"\nFinal Answer: {int(entry_size)} bytes")

solve_lsm_entry_size()
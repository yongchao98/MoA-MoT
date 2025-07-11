import math

def calculate_lsm_io_rate():
    """
    Calculates the minimum total page I/O rate for a leveled LSM tree.
    """
    # 1. Define and standardize inputs
    num_levels = 6  # Number of levels on disk (L0 to L5)
    size_largest_gb = 1  # in Gigabytes
    size_mem_kb = 1      # in Kilobytes
    insert_rate_bps = 16000  # in bytes per second
    page_size_bytes = 2500   # in bytes

    # Convert sizes to bytes
    # The smallest disk level (L0) has the size of the memory buffer flush.
    size_smallest_bytes = size_mem_kb * 1024
    size_largest_bytes = size_largest_gb * 1024**3

    # 2. Calculate the size ratio (T)
    # The number of multiplicative steps from the smallest level (L0) to the largest (L5) is 5.
    num_steps = num_levels - 1
    size_ratio_base = size_largest_bytes / size_smallest_bytes
    size_ratio_t = size_ratio_base**(1 / num_steps)

    # 3. Model and calculate the total byte I/O rate
    # R_io = R_insert * (1 (initial flush) + 2 * (L-1) * T (compactions))
    total_io_rate_bytes_per_sec = insert_rate_bps * (1 + 2 * num_steps * size_ratio_t)

    # 4. Calculate the total page I/O rate
    page_io_rate = total_io_rate_bytes_per_sec / page_size_bytes
    
    # Print the final equation with values plugged in
    print("Equation for Page I/O Rate:")
    print(f"Page I/O Rate = (Insert_Rate * (1 + 2 * (Num_Levels - 1) * Size_Ratio)) / Page_Size")
    print(f"Page I/O Rate = ({insert_rate_bps} * (1 + 2 * ({num_levels} - 1) * {size_ratio_t:.1f})) / {page_size_bytes}")

    # Print the result
    print("\nResult:")
    print(f"The minimum total page I/O rate is: {page_io_rate:.1f} pages/s")

# Execute the calculation
calculate_lsm_io_rate()
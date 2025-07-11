import math

def calculate_lsm_io_rate():
    """
    Calculates the minimum total page I/O rate for an LSM tree based on given parameters.
    """
    # Step 1: Define initial parameters
    num_levels = 6
    size_largest_level_gb = 1.0
    size_mem_buffer_kb = 1.0
    insert_rate_bps = 16000.0
    page_size_bytes = 2500.0

    # Convert sizes to bytes for consistent units
    size_mem_bytes = size_mem_buffer_kb * 1024
    size_largest_level_bytes = size_largest_level_gb * (1024**3)

    print("--- Given Parameters ---")
    print(f"Number of levels (L): {num_levels}")
    print(f"Memory buffer size: {size_mem_bytes} bytes")
    print(f"Largest level size: {size_largest_level_bytes} bytes")
    print(f"Insert rate (R_insert): {insert_rate_bps} bytes/s")
    print(f"Page size (P): {page_size_bytes} bytes")
    print("-" * 26)

    # Step 2: Calculate the size ratio (T)
    # Formula: Size_L(L) = T^L * Size_mem
    # T = (Size_L(L) / Size_mem)^(1/L)
    size_ratio = (size_largest_level_bytes / size_mem_bytes)**(1 / num_levels)
    print("\n--- Step 1: Calculate Size Ratio (T) ---")
    print(f"Equation: T = ({size_largest_level_bytes} / {size_mem_bytes})^(1 / {num_levels})")
    print(f"Result: T = {size_ratio:.4f}")

    # Step 3: Calculate total I/O rate in bytes/s
    # Formula for leveled compaction I/O amplification: (2T + 1)
    # Total I/O Rate = (2T + 1) * R_insert
    io_amplification_factor = 2 * size_ratio + 1
    total_io_rate_bps = io_amplification_factor * insert_rate_bps
    print("\n--- Step 2: Calculate Total I/O Rate (bytes/s) ---")
    print(f"Equation: Total I/O Rate = (2 * {size_ratio:.4f} + 1) * {insert_rate_bps}")
    print(f"Result: Total I/O Rate = {total_io_rate_bps:.4f} bytes/s")

    # Step 4: Calculate minimum total page I/O rate
    # Formula: Page I/O Rate = Total I/O Rate / Page Size
    min_page_io_rate = total_io_rate_bps / page_size_bytes
    print("\n--- Step 3: Calculate Minimum Total Page I/O Rate ---")
    print(f"Final Equation: {total_io_rate_bps:.4f} bytes/s / {page_size_bytes} bytes/page")
    print(f"Final Result: Minimum Page I/O Rate = {min_page_io_rate:.4f} pages/s")

# Execute the calculation
calculate_lsm_io_rate()
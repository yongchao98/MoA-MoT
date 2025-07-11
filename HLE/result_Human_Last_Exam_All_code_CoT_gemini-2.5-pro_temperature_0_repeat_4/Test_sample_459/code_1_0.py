import math

def calculate_lsm_io_rate():
    """
    Calculates the minimum total page I/O rate for an LSM tree based on given parameters.
    """
    # --- Step 1: Define parameters and establish tree structure ---
    # The problem states 6 levels and a memory buffer. We interpret this as
    # 1 memory level (L0) and 5 disk levels (L1 to L5).
    num_total_levels = 6
    num_disk_levels = num_total_levels - 1

    # Given parameters
    largest_level_size_gb = 1.0
    mem_buffer_size_kb = 1.0
    insert_rate_bps = 16000.0  # bytes per second
    page_size_bytes = 2500.0

    # Convert all sizes to bytes for consistency
    largest_level_size_bytes = largest_level_size_gb * (1024**3)
    mem_buffer_size_bytes = mem_buffer_size_kb * 1024

    # --- Step 2: Calculate the Size Ratio (T) ---
    # The size of each level grows by a factor T.
    # Size(L_disk) = Size(L_mem) * T^(num_disk_levels)
    # T = (Size(L_disk) / Size(L_mem))^(1 / num_disk_levels)
    t_pow_l = largest_level_size_bytes / mem_buffer_size_bytes
    size_ratio = t_pow_l**(1.0 / num_disk_levels)

    # --- Step 3 & 4: Calculate the Total I/O Rate in bytes/second ---
    # The total I/O rate is the sum of the memtable flush and all compactions.
    #
    # 1. Memtable Flush (L0 -> L1): Data is written at the insert rate.
    #    I/O_flush = insert_rate_bps (write)
    #
    # 2. Compactions (Li -> L(i+1)): To move data at insert_rate_bps, we must:
    #    - Read from Li: insert_rate_bps
    #    - Read from L(i+1): size_ratio * insert_rate_bps
    #    - Write to L(i+1): (1 + size_ratio) * insert_rate_bps
    #    Total I/O per compaction = (2 + 2 * size_ratio) * insert_rate_bps
    #
    # There are (num_disk_levels - 1) such compactions (L1->L2, ..., L4->L5).
    #
    # Total I/O = I/O_flush + I/O_compactions
    # Total I/O = insert_rate_bps + (num_disk_levels - 1) * (2 + 2 * size_ratio) * insert_rate_bps
    # Total I/O = insert_rate_bps * (1 + (num_disk_levels - 1) * (2 + 2 * size_ratio))

    total_io_bps = insert_rate_bps * (1 + (num_disk_levels - 1) * (2 + 2 * size_ratio))

    # --- Step 5: Convert to Page I/O Rate ---
    page_io_rate = total_io_bps / page_size_bytes

    # --- Final Output ---
    print("Calculation Steps:")
    print(f"1. Tree Structure: {num_total_levels} total levels -> 1 memory buffer (L0) and {num_disk_levels} disk levels (L1-L5).")
    print(f"2. Size Ratio (T): ({largest_level_size_bytes:.0f} / {mem_buffer_size_bytes:.0f})^(1/{num_disk_levels}) = {size_ratio:.4f}")
    print(f"3. Total I/O Rate (bytes/s): {insert_rate_bps:.0f} * (1 + ({num_disk_levels} - 1) * (2 + 2 * {size_ratio:.4f})) = {total_io_bps:.2f} bytes/s")
    print(f"4. Page I/O Rate: {total_io_bps:.2f} bytes/s / {page_size_bytes:.0f} bytes/page = {page_io_rate:.2f} pages/s")
    print("\nFinal Equation:")
    
    # Using integer values for clarity in the final printed equation
    final_equation = (
        f"Minimum total page I/O rate = "
        f"({int(insert_rate_bps)} B/s * (1 + ({num_disk_levels} - 1) * (2 + 2 * {int(round(size_ratio))}.0))) / {int(page_size_bytes)} B/page"
    )
    print(final_equation)
    print(f"= {page_io_rate:.2f} pages/s")


if __name__ == '__main__':
    calculate_lsm_io_rate()
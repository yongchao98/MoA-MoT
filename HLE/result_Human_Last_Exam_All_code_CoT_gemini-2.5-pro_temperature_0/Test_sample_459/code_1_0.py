def calculate_min_lsm_io_rate():
    """
    Calculates the minimum total page I/O rate for an LSM tree.
    """
    # 1. Define the given parameters
    num_total_levels = 6
    insert_rate_bps = 16000  # bytes per second
    page_size_bytes = 2500   # bytes

    # 2. Determine the number of on-disk levels and compactions
    # The 6 levels are assumed to be L0 (memory) and L1-L5 (disk).
    num_on_disk_levels = num_total_levels - 1
    # Compactions occur between on-disk levels: L1->L2, L2->L3, L3->L4, L4->L5
    num_inter_level_compactions = num_on_disk_levels - 1

    print("### Step-by-step Calculation ###")
    print(f"1. System Parameters:")
    print(f"   - Total Levels: {num_total_levels}")
    print(f"   - Insert Rate: {insert_rate_bps} bytes/s")
    print(f"   - Page Size: {page_size_bytes} bytes")
    print(f"   - On-Disk Levels (L1 to L{num_on_disk_levels}): {num_on_disk_levels}")
    print("-" * 30)

    # 3. Calculate the I/O rate for each stage
    print("2. Modeling Minimum I/O Rate:")
    # The initial flush from memory (L0) to the first disk level (L1) is a write operation.
    # The rate of this I/O is equal to the insert rate.
    flush_io_rate_bps = insert_rate_bps
    print(f"   - I/O rate for memory-to-disk flush (L0->L1): {flush_io_rate_bps} bytes/s (write)")

    # For minimum I/O, an inter-level compaction (L_i -> L_{i+1}) involves:
    # - Reading data from L_i at the insert rate.
    # - Writing data to L_{i+1} at the insert rate.
    # The I/O cost per byte is 1 read + 1 write.
    compaction_io_multiplier = 2
    compaction_io_rate_bps = compaction_io_multiplier * insert_rate_bps
    print(f"   - I/O rate for one inter-disk compaction (e.g., L1->L2): {compaction_io_rate_bps} bytes/s (read+write)")
    print(f"   - Number of inter-disk compactions: {num_inter_level_compactions}")
    print("-" * 30)

    # 4. Calculate the total I/O rate in bytes/s
    # Total I/O = (Flush I/O) + (Total Compaction I/O)
    total_io_rate_bps = flush_io_rate_bps + (num_inter_level_compactions * compaction_io_rate_bps)
    
    # Simplified formula: insert_rate * (1 + 2 * num_compactions)
    total_io_factor = 1 + compaction_io_multiplier * num_inter_level_compactions

    print("3. Calculating Total I/O Rate (bytes/s):")
    print(f"   Total I/O Rate = (Flush Rate) + (Num Compactions * Compaction Rate)")
    print(f"   Total I/O Rate = {flush_io_rate_bps} + ({num_inter_level_compactions} * {compaction_io_rate_bps})")
    print(f"   Total I/O Rate = {total_io_rate_bps} bytes/s")
    print("-" * 30)

    # 5. Convert the total I/O rate to pages/s
    total_page_io_rate = total_io_rate_bps / page_size_bytes

    print("4. Final Calculation for Page I/O Rate:")
    print("   The final equation is: (Insert Rate * (1 + 2 * (On-Disk Levels - 1))) / Page Size")
    print(f"   Page I/O Rate = ({insert_rate_bps} * (1 + {compaction_io_multiplier} * {num_inter_level_compactions})) / {page_size_bytes}")
    print(f"   Page I/O Rate = ({insert_rate_bps} * {total_io_factor}) / {page_size_bytes}")
    print(f"   Page I/O Rate = {total_io_rate_bps} / {page_size_bytes}")
    print(f"   Minimum Total Page I/O Rate = {total_page_io_rate} pages/s")

if __name__ == '__main__':
    calculate_min_lsm_io_rate()
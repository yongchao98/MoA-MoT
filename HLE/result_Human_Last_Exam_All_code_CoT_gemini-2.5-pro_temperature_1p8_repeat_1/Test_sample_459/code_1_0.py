import math

def calculate_lsm_io_rate():
    """
    Calculates the minimum total page I/O rate for an LSM tree.
    """
    # 1. Define the given parameters
    num_total_levels = 6
    num_disk_levels = num_total_levels - 1 # L0 is memory, L1-L5 are disk
    mem_buffer_size_kb = 1  # in KB
    largest_level_size_gb = 1  # in GB
    insert_rate_bps = 16000  # in bytes/s
    page_size_bytes = 2500   # in bytes

    # Convert sizes to bytes for consistency
    mem_buffer_size_bytes = mem_buffer_size_kb * 1024
    largest_level_size_bytes = largest_level_size_gb * 1024**3

    print("--- Step 1: Given Parameters ---")
    print(f"Number of levels = {num_total_levels} ({num_disk_levels} on disk)")
    print(f"Memory Buffer Size (L0) = {mem_buffer_size_kb} KB")
    print(f"Largest Level Size (L{num_disk_levels}) = {largest_level_size_gb} GB")
    print(f"Insert Rate = {insert_rate_bps} bytes/s")
    print(f"Page Size = {page_size_bytes} bytes")
    print("-" * 35)

    # 2. Calculate the size ratio (T)
    # Size(L_largest) = T^(num_disk_levels) * Size(L0)
    # T = (Size(L_largest) / Size(L0)) ^ (1 / num_disk_levels)
    size_ratio_t_float = (largest_level_size_bytes / mem_buffer_size_bytes)**(1 / num_disk_levels)
    # Round to the nearest integer as size ratio is typically an integer
    size_ratio_t = int(round(size_ratio_t_float))

    print("--- Step 2: Calculate Size Ratio (T) ---")
    print(f"Equation: T = (Largest Level Size / Memory Buffer Size) ^ (1 / Number of Disk Levels)")
    print(f"T = ({largest_level_size_bytes} / {mem_buffer_size_bytes}) ^ (1 / {num_disk_levels})")
    print(f"T = {size_ratio_t_float:.2f} \u2248 {size_ratio_t}")
    print("-" * 35)

    # 3. Calculate the total I/O rate in bytes/s
    # Total I/O Rate = Initial Flush Write Rate + Compaction I/O Rate
    # Initial Flush Write Rate = insert_rate_bps
    # Compaction I/O for one merge (Li -> L(i+1)):
    #   Read(Li): insert_rate_bps
    #   Read(L(i+1)): T * insert_rate_bps
    #   Write(L(i+1)): T * insert_rate_bps
    #   Total per merge = (1 + 2T) * insert_rate_bps
    # Number of merges between disk levels = num_disk_levels - 1 = 4
    # Total I/O = insert_rate * (1 (initial write) + (num_disk_levels-1) * (1+2T))
    num_inter_disk_merges = num_disk_levels - 1
    total_io_bytes_per_second = insert_rate_bps * (1 + num_inter_disk_merges * (1 + 2 * size_ratio_t))
    
    print("--- Step 3: Calculate Total I/O Rate (Bytes/s) ---")
    print("Formula: Total I/O = Insert Rate * (1 + (Num Disk Levels - 1) * (1 + 2 * T))")
    print(f"Calculation: Total I/O = {insert_rate_bps} * (1 + {num_inter_disk_merges} * (1 + 2 * {size_ratio_t}))")
    final_io_equation_str_part1 = f"{insert_rate_bps} * (1 + {num_inter_disk_merges} * ({1 + 2 * size_ratio_t}))"
    final_io_equation_str_part2 = f"{insert_rate_bps} * ({1 + num_inter_disk_merges * (1 + 2 * size_ratio_t)})"
    print(f"Calculation: Total I/O = {final_io_equation_str_part1} = {final_io_equation_str_part2}")
    print(f"Total I/O Rate = {int(total_io_bytes_per_second)} bytes/s")
    print("-" * 35)
    
    # 4. Convert to page I/O rate
    total_page_io_rate = total_io_bytes_per_second / page_size_bytes

    print("--- Step 4: Convert to Page I/O Rate ---")
    print("Formula: Page I/O Rate = Total I/O Rate (Bytes/s) / Page Size")
    print(f"Final Equation: Total Page I/O Rate = {int(total_io_bytes_per_second)} / {page_size_bytes}")
    print(f"Total Page I/O Rate = {total_page_io_rate:.2f} pages/s")
    print("-" * 35)
    
    return total_page_io_rate

# Run the calculation and print the final answer
final_answer = calculate_lsm_io_rate()
#<<<851.2>>>
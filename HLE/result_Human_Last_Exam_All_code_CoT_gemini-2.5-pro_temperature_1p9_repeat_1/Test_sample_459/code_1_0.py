import math

def calculate_lsm_io_rate():
    """
    Calculates the minimum total page I/O rate for an LSM tree.
    """
    # Step 1: Define the given parameters
    num_levels = 6
    mem_buffer_size_kb = 1
    largest_level_size_gb = 1
    insert_rate_bps = 16000  # bytes per second
    page_size_bytes = 2500

    # Convert sizes to bytes for consistent units
    mem_buffer_size_bytes = mem_buffer_size_kb * 1024
    largest_level_size_bytes = largest_level_size_gb * 1024**3

    print("--- Given Parameters ---")
    print(f"Number of levels (N): {num_levels}")
    print(f"Memory buffer size: {mem_buffer_size_bytes} bytes")
    print(f"Largest level size: {largest_level_size_bytes} bytes")
    print(f"Insert rate: {insert_rate_bps} bytes/s")
    print(f"Page size: {page_size_bytes} bytes")
    print("-" * 26)

    # Step 2: Calculate the size ratio 'k'
    # Formula: L_max = Mem * k^N  =>  k = (L_max / Mem)^(1/N)
    size_ratio_k = (largest_level_size_bytes / mem_buffer_size_bytes) ** (1 / num_levels)
    
    print("\n--- Step 1: Calculate Size Ratio (k) ---")
    print(f"Equation: k = (Largest Level Size / Memory Buffer Size) ^ (1 / Number of Levels)")
    print(f"k = ({largest_level_size_bytes} / {mem_buffer_size_bytes}) ^ (1 / {num_levels})")
    print(f"k = {size_ratio_k:.4f}")
    print("-" * 40)

    # Step 3: Calculate the total I/O rate in bytes per second
    # Formula for leveled compaction: Total I/O = Insert Rate * (N + (N-1) * k)
    # This accounts for N writes and (N-1)*k reads for each byte inserted.
    total_io_bytes_per_sec = insert_rate_bps * (num_levels + (num_levels - 1) * size_ratio_k)

    print("\n--- Step 2: Calculate Total I/O Rate (bytes/s) ---")
    print("Equation: Total I/O Rate = Insert Rate * (N + (N - 1) * k)")
    print(f"Total I/O Rate = {insert_rate_bps} * ({num_levels} + ({num_levels} - 1) * {size_ratio_k:.4f})")
    print(f"Total I/O Rate = {total_io_bytes_per_sec:.4f} bytes/s")
    print("-" * 52)
    
    # Step 4: Convert the I/O rate to pages per second
    page_io_rate = total_io_bytes_per_sec / page_size_bytes

    print("\n--- Step 3: Calculate Total Page I/O Rate (pages/s) ---")
    print("Equation: Page I/O Rate = Total I/O Rate / Page Size")
    print(f"Page I/O Rate = {total_io_bytes_per_sec:.4f} / {page_size_bytes}")
    print(f"Page I/O Rate = {page_io_rate:.4f} pages/s")
    print("-" * 58)

    return page_io_rate

# Run the calculation and print the final answer in the required format
final_answer = calculate_lsm_io_rate()
# print(f"\nFinal Answer: {final_answer:.1f}")
print(f"\n<<< {final_answer:.2f} >>>")

import math

def calculate_lsm_io_rate():
    """
    Calculates the minimum total page I/O rate for an LSM tree based on given parameters.
    """
    # Step 1: Define initial parameters
    num_levels_total = 6
    mem_buffer_size_kb = 1
    largest_level_size_gb = 1
    insert_rate_bps = 16000
    page_size_bytes = 2500

    # Convert sizes to bytes
    mem_buffer_size_bytes = mem_buffer_size_kb * 1024
    largest_level_size_bytes = largest_level_size_gb * 1024 * 1024 * 1024

    print("--- Step 1: System Parameters ---")
    print(f"Total Levels: {num_levels_total}")
    print(f"Memory Buffer Size: {mem_buffer_size_kb} KB ({mem_buffer_size_bytes} bytes)")
    print(f"Largest Level Size: {largest_level_size_gb} GB ({largest_level_size_bytes} bytes)")
    print(f"Insert Rate: {insert_rate_bps} bytes/s")
    print(f"Page Size: {page_size_bytes} bytes\n")

    # Step 2: Calculate the number of disk levels and the size ratio (T)
    # The number of levels on disk is the total minus the one in memory (L0).
    num_disk_levels = num_levels_total - 1
    
    # The size ratio T is derived from Size(L_max) = T^(num_disk_levels) * Size(L0)
    # T^5 = largest_level_size / mem_buffer_size
    size_ratio_base = largest_level_size_bytes / mem_buffer_size_bytes
    T = round(math.pow(size_ratio_base, 1 / num_disk_levels))

    print("--- Step 2: Calculate Size Ratio (T) ---")
    print(f"Number of disk levels (L1 to L{num_disk_levels}): {num_disk_levels}")
    print(f"The size ratio (T) is calculated from: T^{num_disk_levels} = Largest Level Size / Memory Buffer Size")
    print(f"T^{num_disk_levels} = {largest_level_size_bytes} / {mem_buffer_size_bytes} = {size_ratio_base}")
    print(f"T = ({size_ratio_base})^(1/{num_disk_levels}) = {T}\n")

    # Step 3: Calculate the total I/O amplification factor
    # Total Writes per byte = num_disk_levels (one write per level)
    # Total Reads per byte = (num_disk_levels - 1) * (1 read from source + T reads from destination)
    write_amplification = num_disk_levels
    read_amplification = (num_disk_levels - 1) * (1 + T)
    total_io_amplification = write_amplification + read_amplification
    
    print("--- Step 3: Calculate I/O Amplification ---")
    print("For each byte inserted, we calculate the total bytes read and written to disk.")
    print(f"Write Amplification = {num_disk_levels} (1 write for each disk level)")
    print(f"Read Amplification = (num_disk_levels - 1) * (1 + T) = ({num_disk_levels - 1}) * (1 + {T}) = {read_amplification}")
    print(f"Total I/O Amplification = Write Amplification + Read Amplification = {write_amplification} + {read_amplification} = {total_io_amplification}\n")

    # Step 4: Calculate total I/O rate in Bytes/s
    total_io_rate_bps = insert_rate_bps * total_io_amplification
    
    print("--- Step 4: Calculate Total I/O Rate in Bytes/s ---")
    print(f"Total I/O Rate = Insert Rate * Total I/O Amplification")
    print(f"Total I/O Rate = {insert_rate_bps} bytes/s * {total_io_amplification} = {total_io_rate_bps} bytes/s\n")

    # Step 5: Calculate total page I/O rate
    total_page_io_rate = total_io_rate_bps / page_size_bytes

    print("--- Step 5: Calculate Final Page I/O Rate ---")
    print("The final rate is the total I/O rate in bytes divided by the page size.")
    print(f"Final Equation: ({insert_rate_bps} * ({write_amplification} + ({num_disk_levels - 1}) * (1 + {T}))) / {page_size_bytes}")
    print(f"Minimum total page I/O rate = {total_io_rate_bps} bytes/s / {page_size_bytes} bytes/page = {total_page_io_rate} pages/s")


calculate_lsm_io_rate()
<<<467.2>>>
def calculate_ram_test_time():
    """
    This function calculates the time taken for a March RAW test on a 1Mbit RAM.
    """
    # Parameters
    # Test chosen: March RAW, with a complexity of 14N.
    # N is the number of memory bits.
    complexity = 14
    # Size of the RAM in bits (1Mbit = 1,000,000 bits)
    ram_size_bits = 1_000_000
    # Time for one read/write cycle in nanoseconds
    cycle_time_ns = 5

    # Convert cycle time from nanoseconds to seconds for calculation
    cycle_time_s = cycle_time_ns * 1e-9

    # 1. Calculate the total number of read/write operations
    total_operations = complexity * ram_size_bits

    # 2. Calculate the total test time in seconds
    total_time_seconds = total_operations * cycle_time_s

    # 3. Convert the total time to milliseconds
    total_time_ms = total_time_seconds * 1000

    # Print the equation with all the numbers
    print(f"Test Time (ms) = (Complexity * RAM Size * Cycle Time_ns) / 1,000,000")
    print(f"Test Time (ms) = ({complexity} * {ram_size_bits:,} * {cycle_time_ns}) / 1,000,000")
    print(f"Result: {total_time_ms:.0f} ms")


if __name__ == '__main__':
    calculate_ram_test_time()
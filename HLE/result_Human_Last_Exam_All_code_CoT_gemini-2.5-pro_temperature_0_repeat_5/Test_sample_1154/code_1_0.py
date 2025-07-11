def calculate_ram_test_duration():
    """
    Calculates the duration for the RAM test with the highest fault coverage.

    1.  Selects the test with the highest complexity (March RAW, 12N).
    2.  Calculates the total test time for a 1Mbit RAM with a 5ns cycle time.
    3.  Prints the calculation and the final result in milliseconds.
    """
    # Parameters
    # Number of bits in the RAM (N)
    num_bits = 1_000_000
    # Complexity of the chosen test (March RAW is 12N)
    complexity_factor = 12
    # Cycle time for one read/write operation in nanoseconds (tc)
    cycle_time_ns = 5

    # --- Calculations ---
    # Total number of operations
    total_operations = complexity_factor * num_bits
    # Total time in nanoseconds
    total_time_ns = total_operations * cycle_time_ns
    # Convert total time to milliseconds (1 ms = 1,000,000 ns)
    total_time_ms = total_time_ns / 1_000_000

    # --- Output ---
    print("Selected Test: March RAW (highest fault coverage with 12N complexity)")
    print("Calculation for test duration:")
    # Print the equation with all its numbers
    print(f"{complexity_factor} (complexity) * {num_bits:,} (bits) * {cycle_time_ns} (ns/op) = {total_time_ns:,.0f} ns")
    print(f"Total time in milliseconds: {total_time_ms}")

calculate_ram_test_duration()
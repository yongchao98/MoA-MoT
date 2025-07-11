def calculate_ram_test_duration():
    """
    Selects the RAM test with the highest fault coverage and calculates its duration.
    """
    # Step 1: Select the test with the highest fault coverage.
    # Among the given options, March CL and March RAW have the highest complexity (14N),
    # which corresponds to the highest fault coverage. We choose March CL.
    test_name = "March CL"
    complexity_factor = 14

    # Step 2: Define the parameters for the calculation.
    ram_size_bits = 1_000_000  # 1Mbit RAM
    cycle_time_ns = 5          # 5 nanoseconds per read/write cycle
    ns_in_a_ms = 1_000_000     # Conversion factor from nanoseconds to milliseconds

    # Step 3: Calculate the total test duration.
    # Total time = (Complexity * RAM Size * Cycle Time)
    total_time_ns = complexity_factor * ram_size_bits * cycle_time_ns
    total_time_ms = total_time_ns / ns_in_a_ms

    # Step 4: Print the details of the calculation and the final result.
    print(f"Chosen Test: {test_name} (Complexity: {complexity_factor}N)")
    print(f"RAM Size (N): {ram_size_bits:,} bits")
    print(f"Cycle Time (tc): {cycle_time_ns} ns")
    print("\nCalculating the total test duration:")
    # The user requested to output each number in the final equation
    print(f"({complexity_factor} operations/bit * {ram_size_bits:,} bits * {cycle_time_ns} ns/operation) / {ns_in_a_ms:,} ns/ms")
    
    print(f"\nResult: {total_time_ms} milliseconds")


if __name__ == "__main__":
    calculate_ram_test_duration()

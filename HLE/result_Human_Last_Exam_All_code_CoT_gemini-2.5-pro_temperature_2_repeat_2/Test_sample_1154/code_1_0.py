def calculate_ram_test_time():
    """
    Calculates the duration of a March RAW RAM test.

    This function determines the RAM test with the highest fault coverage from a
    predefined list, and then calculates the time it would take to run this
    test on a 1Mbit RAM with a 5ns cycle time.
    """

    # --- Step 1: Define parameters ---
    # The test with the highest fault coverage from the list is March RAW.
    # The complexity of the March RAW test is 13N, where N is the number of bits.
    test_name = "March RAW"
    complexity_factor = 13

    # RAM size in bits
    ram_size_bits = 1000000

    # Time for one read/write cycle in nanoseconds
    cycle_time_ns = 5

    # Conversion factor from nanoseconds to milliseconds
    ns_to_ms_conversion = 1000000

    # --- Step 2: Calculate the total test time ---
    total_operations = complexity_factor * ram_size_bits
    total_time_ns = total_operations * cycle_time_ns
    total_time_ms = total_time_ns / ns_to_ms_conversion

    # --- Step 3: Print the results ---
    print(f"Selected test with the highest fault coverage: {test_name} (Complexity: {complexity_factor}N)")
    print(f"RAM Size (N): {ram_size_bits} bits")
    print(f"Cycle Time (tc): {cycle_time_ns} ns")
    print("\nCalculating the total test duration in milliseconds:")
    # The user requested to output each number in the final equation.
    print(f"Formula: (Complexity Factor * RAM Size * Cycle Time) / ns-to-ms-conversion")
    print(f"Calculation: ({complexity_factor} * {ram_size_bits} * {cycle_time_ns}) / {ns_to_ms_conversion}")

    # Final result
    print(f"\nResult: {total_time_ms} ms")


if __name__ == "__main__":
    calculate_ram_test_time()
    # The final number representing the duration in milliseconds
    print("<<<65.0>>>")

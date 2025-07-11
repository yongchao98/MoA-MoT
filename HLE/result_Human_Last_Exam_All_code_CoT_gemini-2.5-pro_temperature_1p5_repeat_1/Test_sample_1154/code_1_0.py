def calculate_ram_test_time():
    """
    Calculates the duration of a March RAW test on a 1Mbit RAM.
    """
    # Step 1: Define the parameters for the calculation.
    # The chosen test is March RAW, which has the highest fault coverage.
    # Its complexity is 14N, meaning 14 operations per bit.
    test_name = "March RAW"
    complexity_factor = 14

    # RAM size is 1 Mbit (1,000,000 bits).
    ram_size_bits = 1_000_000

    # Cycle time for one read or write operation is 5 nanoseconds.
    cycle_time_ns = 5

    # Step 2: Calculate the total test time.
    # Total operations = RAM size * Complexity factor
    total_operations = ram_size_bits * complexity_factor

    # Total time in nanoseconds = Total operations * Cycle time
    total_time_ns = total_operations * cycle_time_ns

    # Convert the time from nanoseconds to milliseconds (1 ms = 1,000,000 ns).
    total_time_ms = total_time_ns / 1_000_000

    # Step 3: Print the results, including the numbers used in the final equation.
    print(f"Chosen Test: {test_name} (Complexity: {complexity_factor}N)")
    print("\nCalculation of test duration:")
    print(f"Total time = (RAM size) * (Complexity factor) * (Cycle time)")
    # The final equation with all numbers substituted
    print(f"{ram_size_bits} * {complexity_factor} * {cycle_time_ns} ns = {total_time_ns} ns")
    print(f"{total_time_ns} ns is equal to {int(total_time_ms)} milliseconds.")
    print(f"\nThe final answer is: {int(total_time_ms)}")

calculate_ram_test_time()
def calculate_ram_test_time():
    """
    Calculates the time for a March RAW test on a 1Mbit RAM.
    """
    # The chosen test is March RAW, which has the highest complexity (13N)
    # and therefore the highest fault coverage.
    complexity_factor = 13

    # Number of bits in the RAM (1 Mbit)
    n_bits = 1_000_000

    # Time for one read/write cycle in nanoseconds
    cycle_time_ns = 5

    # Conversion factor from nanoseconds to milliseconds
    ns_in_one_ms = 1_000_000

    # Calculate the total test duration in milliseconds
    test_duration_ms = (complexity_factor * n_bits * cycle_time_ns) / ns_in_one_ms

    # Print the explanation and the final equation as requested
    print("Test with highest fault coverage: March RAW (13N)")
    print("Calculation for a 1Mbit RAM with a 5ns cycle time:")
    print(f"{complexity_factor} * {n_bits} * {cycle_time_ns} / {ns_in_one_ms} = {test_duration_ms}")

# Execute the function
calculate_ram_test_time()
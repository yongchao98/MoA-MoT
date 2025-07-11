def calculate_ram_test_time():
    """
    Calculates the time for a March RAW test on a 1Mbit RAM.
    """
    # 1. Define the parameters
    # RAM size in bits
    ram_size_bits = 1000000
    # Time for one read/write cycle in nanoseconds
    cycle_time_ns = 5
    # The test with the highest fault coverage is March RAW, with a complexity of 14N.
    # This means 14 operations are performed per bit.
    test_complexity_factor = 14

    # 2. Calculate the total time in nanoseconds
    total_time_ns = ram_size_bits * test_complexity_factor * cycle_time_ns

    # 3. Convert the total time to milliseconds (1 ms = 1,000,000 ns)
    total_time_ms = total_time_ns / 1000000

    # 4. Print the explanation and the result
    print("Step 1: The test with the highest fault coverage is March RAW (14N).")
    print("Step 2: Calculate the test duration for a 1Mbit RAM with a 5ns cycle time.")
    print("\nEquation:")
    print(f"Test Duration (ms) = (RAM Size * Complexity Factor * Cycle Time) / 1,000,000")
    # Output each number in the final equation
    print(f"Test Duration (ms) = ({ram_size_bits} * {test_complexity_factor} * {cycle_time_ns}) / 1000000")
    print(f"\nResult: {total_time_ms} ms")

calculate_ram_test_time()
def calculate_ram_test_time():
    """
    Calculates the time taken for the March RAW test on a 1Mbit RAM.
    """
    # 1. Define the parameters for the calculation.

    # The chosen test is March RAW, which has the highest complexity (13N)
    # among the given options, indicating the highest fault coverage.
    complexity = 13  # operations per bit for March RAW

    # The RAM size is 1Mbit.
    ram_size_bits = 1000000  # 1 million bits

    # The time for one read/write cycle.
    cycle_time_ns = 5  # 5 nanoseconds

    # Conversion factor: there are 1,000,000 nanoseconds in 1 millisecond.
    ns_to_ms_divisor = 1000000

    # 2. Perform the calculation.
    total_time_ms = (complexity * ram_size_bits * cycle_time_ns) / ns_to_ms_divisor

    # 3. Print the results, including the equation with all its numbers.
    print("Test with highest fault coverage: March RAW (13N complexity)")
    print("-" * 50)
    print("Calculation for a 1Mbit RAM with a 5ns cycle time:")
    
    # Show the equation with all the numbers, as requested.
    print(f"({complexity} operations/bit * {ram_size_bits:,} bits * {cycle_time_ns} ns/operation) / {ns_to_ms_divisor:,} ns/ms")

    # Print the final result.
    print(f"= {total_time_ms} ms")

calculate_ram_test_time()
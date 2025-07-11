def calculate_ram_test_duration():
    """
    Calculates the duration for a RAM test with the highest fault coverage from a given list.
    """
    # Step 1: Identify the test with the highest fault coverage.
    # Among the options, March CL and March RAW have the highest complexity (12N),
    # indicating the highest fault coverage. We will choose March CL.
    test_name = "March CL"
    complexity = 12

    # Step 2: Define the given parameters.
    num_bits = 1000000      # 1 Mbit = 1,000,000 bits
    cycle_time_ns = 5       # 5 ns per read/write cycle

    # Step 3: Calculate the total number of operations.
    total_operations = complexity * num_bits

    # Step 4: Calculate the total time in nanoseconds.
    total_time_ns = total_operations * cycle_time_ns

    # Step 5: Convert the total time to milliseconds.
    # 1 ms = 1,000,000 ns
    conversion_factor = 1000000
    total_time_ms = total_time_ns / conversion_factor

    # Step 6: Print the detailed calculation process.
    print(f"The RAM test with the highest fault coverage from the list is {test_name}, with a complexity of {complexity}N.")
    print("\nCalculation of the test duration:")
    print(f"Total time = (Complexity * Number of bits) * Time per operation")
    print(f"Total time = ({complexity} * {num_bits}) * {cycle_time_ns} ns")
    print(f"Total time = {total_time_ns} ns")
    print("\nConverting nanoseconds to milliseconds:")
    print(f"Total time = {total_time_ns} ns / {conversion_factor} ns/ms")
    print(f"Duration = {total_time_ms} ms")

calculate_ram_test_duration()

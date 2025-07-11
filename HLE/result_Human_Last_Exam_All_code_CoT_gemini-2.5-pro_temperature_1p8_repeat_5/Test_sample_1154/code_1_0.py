def calculate_ram_test_duration():
    """
    Calculates the duration of a March RAW test on a 1Mbit RAM.
    """
    # Step 1: Define the parameters for the calculation.
    # The chosen test is March RAW, which has the highest fault coverage.
    # Its complexity is 14N, meaning 14 operations per bit.
    test_name = "March RAW"
    complexity_factor = 14
    
    # RAM size in bits (1Mbit)
    ram_size_bits = 1_000_000
    
    # Time for one read/write cycle in nanoseconds
    time_per_cycle_ns = 5

    # Step 2: Perform the calculation.
    # Calculate the total number of operations.
    total_operations = complexity_factor * ram_size_bits

    # Calculate the total test time in nanoseconds.
    total_time_ns = total_operations * time_per_cycle_ns

    # Convert the total time to milliseconds (1 ms = 1,000,000 ns).
    total_time_ms = total_time_ns / 1_000_000
    
    # Step 3: Print the results following the requested format.
    print(f"Chosen Test with Highest Fault Coverage: {test_name}")
    print(f"Test Complexity: {complexity_factor}N")
    print(f"RAM size (N): {ram_size_bits:,} bits")
    print(f"Time per cycle (tc): {time_per_cycle_ns} ns")
    print("\nCalculation Steps:")
    print(f"Total Operations = {complexity_factor} * {ram_size_bits:,} = {total_operations:,}")
    print(f"Total Time (ns) = {total_operations:,} operations * {time_per_cycle_ns} ns/operation = {total_time_ns:,} ns")
    print(f"Total Time (ms) = {total_time_ns:,} ns / 1,000,000 ns/ms = {int(total_time_ms)} ms")
    print("\n---")
    print(f"The total time taken for the test is {int(total_time_ms)} milliseconds.")
    
# Execute the function
calculate_ram_test_duration()
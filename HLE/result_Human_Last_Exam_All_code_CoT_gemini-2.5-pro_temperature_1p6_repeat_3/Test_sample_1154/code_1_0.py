def calculate_ram_test_duration():
    """
    Calculates the duration of the March RAW test for a 1Mbit RAM.
    """
    # Parameters
    # The test with the highest fault coverage from the list is March RAW.
    # Its complexity is 16N, where N is the number of bits.
    test_name = "March RAW"
    complexity_factor = 16
    
    # RAM size in bits (1 Mbit)
    ram_size_bits = 1_000_000
    
    # Cycle time for one read/write operation in nanoseconds
    cycle_time_ns = 5

    # --- Calculation ---
    
    # 1. Total number of operations (reads/writes)
    total_operations = complexity_factor * ram_size_bits
    
    # 2. Total time in nanoseconds
    total_time_ns = total_operations * cycle_time_ns
    
    # 3. Total time in seconds (1 second = 1,000,000,000 nanoseconds)
    total_time_s = total_time_ns / 1_000_000_000
    
    # 4. Total time in milliseconds (1 second = 1,000 milliseconds)
    total_time_ms = total_time_s * 1000

    # --- Output ---
    
    print(f"The selected test with the highest fault coverage is {test_name}.\n")
    print("The calculation for the test duration is as follows:\n")
    
    print("Formula: Total Time = (RAM Size) * (Complexity Factor) * (Cycle Time)\n")
    
    # Output each number in the final equation
    print("Equation with values:")
    print(f"Total Time (ms) = {ram_size_bits:,} bits * {complexity_factor} * {cycle_time_ns} ns\n")

    print(f"1. Total operations = {complexity_factor} * {ram_size_bits:,} = {total_operations:,} operations")
    print(f"2. Total time in seconds = {total_operations:,} * {cycle_time_ns} ns = {total_time_s} s")
    print(f"3. Total time in milliseconds = {total_time_s} s * 1000 = {int(total_time_ms)} ms\n")

    print("The final result for the duration of the test in milliseconds is:")
    print(int(total_time_ms))

calculate_ram_test_duration()
def calculate_ram_test_time():
    """
    This script calculates the time taken for a March RAW test on a 1Mbit RAM.
    """
    # Step 1: Define the parameters for the calculation.
    # The test with the highest fault coverage from the list is March RAW.
    # Its complexity is 13N, meaning 13 operations per bit.
    test_complexity = 13
    
    # RAM size is 1 Mbit, which is 1,000,000 bits.
    ram_size_bits = 1_000_000
    
    # One read/write cycle takes 5 nanoseconds.
    time_per_cycle_ns = 5
    
    # Step 2: Calculate the total number of operations.
    total_operations = ram_size_bits * test_complexity
    
    # Step 3: Calculate the total time in nanoseconds.
    total_time_ns = total_operations * time_per_cycle_ns
    
    # Step 4: Convert the total time to milliseconds.
    # 1 millisecond = 1,000,000 nanoseconds.
    total_time_ms = total_time_ns / 1_000_000
    
    # Print the explanation and the breakdown of the calculation.
    print("Chosen Test: March RAW (highest fault coverage with 13N complexity)")
    print("\n--- Calculation Breakdown ---")
    print(f"Total Operations = RAM Size (bits) * Complexity")
    print(f"                 = {ram_size_bits:,} * {test_complexity} = {total_operations:,} operations")
    
    print(f"\nTotal Time (ns)  = Total Operations * Time per Cycle (ns)")
    print(f"                 = {total_operations:,} * {time_per_cycle_ns} = {total_time_ns:,} ns")

    print(f"\nTotal Time (ms)  = Total Time (ns) / 1,000,000")
    print(f"                 = {total_time_ns:,} / 1,000,000 = {total_time_ms} ms")
    
    # The final question asks to return only the number.
    print("\nFinal Answer (milliseconds):")
    print(total_time_ms)

# Execute the function to perform the calculation and print the results.
calculate_ram_test_time()
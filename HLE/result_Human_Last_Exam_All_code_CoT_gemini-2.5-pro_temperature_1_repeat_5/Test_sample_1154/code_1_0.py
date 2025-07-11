def calculate_ram_test_time():
    """
    This function calculates the time taken for a March RAW test on a 1Mbit RAM.
    """
    # Step 1: Define the parameters for the calculation.
    # The chosen test is March RAW, which has the highest fault coverage among the list.
    # Its complexity is 14N, meaning 14 operations per bit.
    complexity_factor = 14
    
    # RAM size is 1Mbit (1,000,000 bits).
    ram_size_bits = 1000000
    
    # Time for one read/write cycle is 5 nanoseconds.
    cycle_time_ns = 5
    
    # Step 2: Calculate the total number of operations.
    # Total Operations = Complexity * RAM Size
    total_operations = complexity_factor * ram_size_bits
    
    # Step 3: Calculate the total time in nanoseconds.
    # Total Time (ns) = Total Operations * Cycle Time (ns)
    total_time_ns = total_operations * cycle_time_ns
    
    # Step 4: Convert the total time to milliseconds.
    # 1 millisecond = 1,000,000 nanoseconds
    total_time_ms = total_time_ns / 1000000
    
    # Print the explanation with the numbers used in the final equation.
    print("Chosen Test: March RAW (Complexity: 14N)")
    print("\n--- Calculation ---")
    print(f"Total Time = (Complexity Factor * RAM Size) * Cycle Time")
    print(f"Equation: ({complexity_factor} * {ram_size_bits:,}) * {cycle_time_ns} ns")
    print(f"Result in ns: {total_time_ns:,} ns")
    print(f"Result in ms: {total_time_ns:,} ns / 1,000,000 = {total_time_ms} ms")
    
    # Return the final numerical answer.
    return total_time_ms

# Execute the function to print the result.
# The final result is the duration in milliseconds.
final_answer = calculate_ram_test_time()

# The problem asks for a single number as the duration.
# print(f"\nFinal Answer (duration in ms): {final_answer}")
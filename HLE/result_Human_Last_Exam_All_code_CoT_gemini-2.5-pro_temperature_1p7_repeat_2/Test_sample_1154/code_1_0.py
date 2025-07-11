import math

def calculate_ram_test_time():
    """
    This function calculates the time for a RAM test based on its complexity.

    The steps are:
    1.  Choose the test with the highest fault coverage. Based on standard March test algorithms,
        March CL is chosen for its high complexity (16N), which correlates with high fault coverage.
    2.  Define the given parameters.
    3.  Calculate the total test time.
    4.  Print the breakdown of the calculation and the final result in milliseconds.
    """
    
    # Parameters
    # From the list, March CL is selected as it has one of the highest complexities
    # and therefore highest fault coverage. Its complexity is typically 16N.
    complexity_factor = 16
    # RAM size in bits (1 Mbit)
    ram_size_n = 1_000_000
    # Time per read/write cycle in nanoseconds
    time_per_cycle_ns = 5

    # Calculation
    # Total operations = complexity_factor * N
    total_operations = complexity_factor * ram_size_n
    
    # Total time in nanoseconds = total_operations * time_per_cycle
    total_time_ns = total_operations * time_per_cycle_ns
    
    # Convert total time from nanoseconds to milliseconds (1 ms = 1,000,000 ns)
    total_time_ms = total_time_ns / 1_000_000

    # Print the explanation and the final equation with values
    print("--- RAM Test Time Calculation ---")
    print(f"Selected Test: March CL (highest fault coverage)")
    print(f"Algorithm Complexity: {complexity_factor}N")
    print("-" * 20)
    print("Equation: Time = Complexity_Factor * RAM_Size * Time_per_Cycle")
    print(f"Values: Time = {complexity_factor} * {ram_size_n:,} bits * {time_per_cycle_ns} ns")
    print("-" * 20)
    
    # Print the result
    print(f"Total calculated time: {total_time_ns:,} ns")
    print(f"Total time in milliseconds: {total_time_ms}")

# Run the calculation
calculate_ram_test_time()
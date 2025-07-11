import sys

def solve():
    """
    This script identifies the RAM test with the highest fault coverage and calculates its duration.
    
    Plan:
    1. Choose the RAM test with the highest fault coverage from the provided list.
       March CL (12N complexity) is selected for its superior fault coverage over other
       listed algorithms like March C- (10N) or March X (6N).
    2. Define the parameters for the calculation:
       - Test complexity factor (k): 12 for March CL.
       - RAM size (N): 1 Mbit = 1,000,000 bits.
       - Cycle time (tc): 5 nanoseconds (ns).
    3. Calculate the total test duration in milliseconds.
    4. Print the output showing the chosen test, the full equation, and the final result.
    """
    # Step 1: Explanation
    print("Chosen Test: March CL")
    print("Reason: It has a complexity of 12N and offers the highest fault coverage among the listed options, detecting a wide array of fault types.\n")
    
    # Step 2: Define parameters
    # N = number of bits in the RAM
    N_bits = 1_000_000
    # k = complexity factor for the chosen test (12N for March CL)
    k_complexity = 12
    # tc_ns = time for one read/write cycle in nanoseconds
    tc_ns = 5
    
    # Step 3: Calculate the total test duration
    total_operations = k_complexity * N_bits
    total_time_ns = total_operations * tc_ns
    # Convert nanoseconds to milliseconds (1 ms = 1,000,000 ns)
    total_time_ms = total_time_ns / 1_000_000
    
    # Step 4: Print the final calculation and result
    print("Calculation of the test duration:")
    # Print the full equation with all numbers
    print(f"Test Duration = {k_complexity} * {N_bits:,} bits * {tc_ns} ns")
    print(f"Total time = {total_time_ns:,} ns")
    print(f"Result: The total time for the test is {total_time_ms} milliseconds.")

solve()
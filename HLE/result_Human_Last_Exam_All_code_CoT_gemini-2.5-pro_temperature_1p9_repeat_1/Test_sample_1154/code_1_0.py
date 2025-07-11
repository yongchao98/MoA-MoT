def calculate_ram_test_time():
    """
    This script selects the RAM test with the highest fault coverage from the list
    and calculates its execution time on a 1Mbit RAM.
    """
    
    # 1. Select the test and define parameters
    test_name = "March RAW"
    fault_coverage_reason = "It covers the same faults as March C- plus Read Destructive Faults (RDFs)."
    
    N = 1_000_000  # RAM size in bits (1Mbit)
    k = 14         # Algorithmic complexity for March RAW (14N)
    tc_ns = 5      # Cycle time in nanoseconds

    print(f"Selected Test: {test_name}")
    print(f"Reason for Selection: Highest fault coverage. {fault_coverage_reason}")
    print("-" * 30)
    
    # 2. Perform the calculation
    total_operations = N * k
    total_time_ns = total_operations * tc_ns
    total_time_ms = total_time_ns / 1_000_000 # Convert nanoseconds to milliseconds

    # 3. Print the results including the equation
    print("Calculation Steps:")
    print(f"Test Time = RAM Size (N) * Complexity (k) * Cycle Time (tc)")
    print("\nThe final equation is:")
    # Printing each number in the final equation as requested
    print(f"Time = {N:,} bits * {k} operations/bit * {tc_ns} ns/operation")
    print(f"Time = {total_operations:,} total operations * {tc_ns} ns/operation")
    print(f"Time = {total_time_ns:,} ns")
    print("-" * 30)
    print(f"Total test duration in milliseconds: {int(total_time_ms)}")

calculate_ram_test_time()
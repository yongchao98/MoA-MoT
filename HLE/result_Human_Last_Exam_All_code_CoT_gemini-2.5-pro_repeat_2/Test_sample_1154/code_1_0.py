def calculate_ram_test_time():
    """
    This script identifies the RAM test with the highest fault coverage from the list
    and calculates its execution time for a 1Mbit RAM.
    """
    # Test chosen: March RAW, for its highest fault coverage (complexity 14N).
    test_name = "March RAW"
    complexity = 14

    # Given parameters
    ram_size_bits = 1_000_000
    cycle_time_ns = 5

    # --- Calculation ---
    # 1. Total operations = complexity * number of bits
    total_operations = complexity * ram_size_bits

    # 2. Total time in nanoseconds = total operations * time per operation
    total_time_ns = total_operations * cycle_time_ns

    # 3. Convert nanoseconds to milliseconds (1 ms = 1,000,000 ns)
    total_time_ms = total_time_ns / 1_000_000

    # --- Output ---
    print(f"The chosen test with the highest fault coverage is: {test_name}")
    print(f"This test has a complexity of {complexity}N.")
    print("\n--- Calculation Details ---")
    print(f"RAM Size (N):         {ram_size_bits:,} bits")
    print(f"Complexity Factor:    {complexity}")
    print(f"Cycle Time (tc):      {cycle_time_ns} ns")
    print("\n--- Final Equation ---")
    print("Test Time (ms) = (RAM Size * Complexity Factor * Cycle Time) / 1,000,000")
    print(f"Test Time (ms) = ({ram_size_bits:,} * {complexity} * {cycle_time_ns}) / 1,000,000")
    print(f"\nResult: The total test duration is {total_time_ms} milliseconds.")

if __name__ == "__main__":
    calculate_ram_test_time()
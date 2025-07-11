def calculate_ram_test_duration():
    """
    Calculates the duration of a March RAW test for a 1Mbit RAM.
    """
    # Test with the highest fault coverage from the list is March RAW.
    # Its complexity is 14N.
    complexity_factor = 14

    # Number of bits in a 1Mbit RAM
    num_bits = 1000000

    # Time per read/write cycle in nanoseconds
    cycle_time_ns = 5

    # --- Calculation ---
    # Total time in nanoseconds
    total_time_ns = complexity_factor * num_bits * cycle_time_ns

    # Convert total time to milliseconds
    # 1 ms = 1,000,000 ns
    total_time_ms = total_time_ns / 1000000.0

    # --- Output ---
    print("Chosen Test: March RAW (highest fault coverage from the list)")
    print(f"Test Complexity: {complexity_factor}N")
    print("\nCalculating test duration for a 1Mbit RAM with a 5ns cycle time:")
    
    # As requested, printing each number in the final equation
    print("\nEquation: Test Duration = Complexity Factor * Number of Bits * Cycle Time")
    print(f"Values:   Test Duration = {complexity_factor} * {num_bits:,} * {cycle_time_ns} ns")
    
    print(f"\nTotal duration in nanoseconds: {int(total_time_ns):,}")
    print(f"Total duration in milliseconds: {total_time_ms}")

if __name__ == '__main__':
    calculate_ram_test_duration()
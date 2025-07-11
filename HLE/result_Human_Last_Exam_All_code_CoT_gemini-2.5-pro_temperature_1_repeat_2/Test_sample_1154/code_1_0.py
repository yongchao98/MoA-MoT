def calculate_ram_test_time():
    """
    Identifies the RAM test with the highest fault coverage from a given list
    and calculates the time it would take to run on a 1Mbit RAM.
    """
    # Step 1: Identify the test with the highest fault coverage.
    # The complexity of a memory test (k in k*N operations) is a good
    # indicator of its fault coverage.
    #
    # Test Complexities:
    # - MSCAN: 4N
    # - MATS: 4N
    # - March X: 6N
    # - March Y: 8N
    # - March C-: 10N
    # - March C: 11N
    # - March RAW: 12N
    # - March CL: 12N
    #
    # March RAW and March CL have the highest complexity (12N), and thus the highest fault coverage.
    # We will choose March RAW for this calculation.
    test_name = "March RAW"
    complexity_factor = 12  # for a 12N test

    # Step 2: Define parameters for the calculation.
    ram_size_bits = 1_000_000  # 1 Mbit RAM
    cycle_time_ns = 5         # 5 nanoseconds per read/write cycle
    ns_in_one_ms = 1_000_000  # Conversion factor from nanoseconds to milliseconds

    # Step 3: Calculate the total test time.
    total_time_ns = complexity_factor * ram_size_bits * cycle_time_ns
    total_time_ms = total_time_ns / ns_in_one_ms

    # Step 4: Print the explanation, equation, and result.
    print(f"The RAM test chosen for its highest fault coverage is: {test_name} (Complexity: {complexity_factor}N)")
    print("\nThe calculation for the test duration is as follows:")
    print("Test Duration (ms) = (Complexity * RAM Size in bits * Cycle Time in ns) / (ns per ms)")
    
    # Print the equation with the actual numbers
    print("\nFinal Equation:")
    print(f"Test Duration (ms) = ({complexity_factor} * {ram_size_bits} * {cycle_time_ns}) / {ns_in_one_ms}")

    print(f"\nCalculated Test Duration: {total_time_ms} ms")

if __name__ == '__main__':
    calculate_ram_test_time()

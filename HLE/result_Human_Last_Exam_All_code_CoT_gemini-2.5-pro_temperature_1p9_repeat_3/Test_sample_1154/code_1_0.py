def calculate_ram_test_duration():
    """
    Calculates the duration for the RAM test with the highest fault coverage.
    """
    # Step 1: Choose the test with the highest fault coverage.
    # From the list (MSCAN, MATS, March X, Y, RAW, C, C-, CL), March RAW
    # has the highest complexity (13N) and thus the highest fault coverage,
    # as it detects complex faults like Read Destructive Faults.
    test_name = "March RAW"
    complexity_factor = 13

    # Step 2: Define the given parameters.
    # RAM size in bits (N)
    N = 1_000_000
    # Read/write cycle time in nanoseconds (tc)
    tc_ns = 5
    # Conversion factor from nanoseconds to milliseconds
    ns_to_ms_conversion = 1_000_000

    # Step 3: Calculate the total test duration.
    # Total time = (Complexity Factor * RAM Size * Cycle Time)
    total_time_ns = complexity_factor * N * tc_ns
    total_time_ms = total_time_ns / ns_to_ms_conversion

    print(f"Selected Test: {test_name} (Complexity: {complexity_factor}N)")
    print(f"RAM Size (N): {N:,} bits")
    print(f"Cycle Time (tc): {tc_ns} ns")
    print("-" * 30)
    print("The total test duration is calculated as:")
    print(f"({complexity_factor} operations/bit * {N:,} bits * {tc_ns} ns/op) / {ns_to_ms_conversion:,} ns/ms")
    print(f"Result: {total_time_ms} ms")

calculate_ram_test_duration()
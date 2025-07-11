def calculate_ram_test_time():
    """
    Calculates the time taken for a March RAW test on a 1Mbit RAM.
    """
    # Step 1: Define the parameters
    # The test with the highest fault coverage is March RAW, with 14N complexity.
    complexity_factor = 14
    # N is the number of bits in the RAM
    num_bits = 1_000_000
    # tc is the time for one read/write cycle in nanoseconds
    time_per_cycle_ns = 5
    # Conversion factor from nanoseconds to milliseconds
    ns_to_ms_conversion = 1_000_000

    print("Chosen Test: March RAW (Complexity: 14N) for its high fault coverage.")
    print("-" * 30)
    
    # Step 2 & 3: Calculate the total time in nanoseconds
    total_time_ns = complexity_factor * num_bits * time_per_cycle_ns

    # Step 4: Convert the time to milliseconds
    total_time_ms = total_time_ns / ns_to_ms_conversion

    # Output the final equation with all numbers
    print("Calculation Equation:")
    print(f"Time (ms) = ({complexity_factor} operations/bit * {num_bits:,} bits * {time_per_cycle_ns} ns/operation) / {ns_to_ms_conversion:,} ns/ms")
    print("-" * 30)

    # Print the final result
    print(f"Total Test Duration: {total_time_ms} ms")

if __name__ == "__main__":
    calculate_ram_test_time()
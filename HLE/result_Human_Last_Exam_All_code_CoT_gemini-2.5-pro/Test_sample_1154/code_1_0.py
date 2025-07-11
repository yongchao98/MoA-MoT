def calculate_ram_test_duration():
    """
    Calculates the duration of the March RAW test for a 1Mbit RAM.
    """
    # Step 1: Define test parameters
    # The March RAW test is chosen for its high fault coverage (13N complexity).
    complexity_factor = 13
    ram_size_bits = 1_000_000
    cycle_time_ns = 5

    # Step 2: Calculate the total number of operations
    total_operations = complexity_factor * ram_size_bits

    # Step 3: Calculate the total time in nanoseconds
    total_time_ns = total_operations * cycle_time_ns

    # Step 4: Convert the total time to milliseconds
    # 1 millisecond = 1,000,000 nanoseconds
    total_time_ms = total_time_ns / 1_000_000

    # Step 5: Print the breakdown of the calculation
    print(f"Chosen Test: March RAW (Complexity: {complexity_factor}N)")
    print("Calculation:")
    print(f"Total Time = (Test Complexity * RAM Size * Cycle Time)")
    print(f"Total Time (ms) = ({complexity_factor} * {ram_size_bits:,} * {cycle_time_ns}) / 1,000,000")
    print(f"Total Time (ms) = ({total_operations:,} * {cycle_time_ns}) / 1,000,000")
    print(f"Total Time (ms) = {total_time_ns:,} / 1,000,000")
    print(f"Final Duration: {int(total_time_ms)} milliseconds")


if __name__ == "__main__":
    calculate_ram_test_duration()
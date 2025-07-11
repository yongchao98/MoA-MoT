def calculate_ram_test_duration():
    """
    Calculates the duration of the March CL RAM test for a 1Mbit RAM.
    """
    # Step 1: Define the parameters based on the problem statement.
    # The test with the highest fault coverage is March CL, with a complexity of 12N.
    complexity_factor = 12
    ram_size_in_bits = 1_000_000
    time_per_cycle_ns = 5
    
    # Conversion factor from nanoseconds to milliseconds.
    ns_in_one_ms = 1_000_000

    # Step 2: Calculate the total time.
    # Total time in ns = Complexity * RAM size * Time per cycle
    total_time_in_ns = complexity_factor * ram_size_in_bits * time_per_cycle_ns
    
    # Convert total time to milliseconds.
    total_time_in_ms = total_time_in_ns / ns_in_one_ms

    # Step 3: Print the explanation and the final equation with numbers.
    print("The selected test with the highest fault coverage is March CL (12N complexity).")
    print("The total test duration is calculated using the following formula:")
    print("Time (ms) = (Complexity Factor * RAM Size (bits) * Time per Cycle (ns)) / Nanoseconds_per_Millisecond\n")
    print("Plugging in the numbers:")
    # Using f-string formatting to display the equation with its values
    print(f"Time (ms) = ({complexity_factor} * {ram_size_in_bits:,} * {time_per_cycle_ns}) / {ns_in_one_ms:,}")
    print(f"Result = {total_time_in_ms}")

# Execute the function to see the output.
calculate_ram_test_duration()

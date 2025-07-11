def calculate_ram_test_time():
    """
    Calculates the time for a RAM test based on its complexity.
    """
    # 1. Define problem parameters.
    ram_size_bits = 1_000_000
    time_per_cycle_ns = 5.0

    # 2. Choose the test with the highest fault coverage.
    # March RAW is chosen due to its high complexity (12N), which provides
    # broader fault coverage than other listed tests like March C- (10N) or March Y (8N).
    test_name = "March RAW"
    test_complexity = 12

    # 3. Calculate the total test time.
    total_operations = ram_size_bits * test_complexity
    total_time_ns = total_operations * time_per_cycle_ns

    # 4. Convert the result to milliseconds for the final answer.
    # 1 millisecond = 1,000,000 nanoseconds
    total_time_ms = total_time_ns / 1_000_000

    # 5. Print the explanation and step-by-step equation.
    print(f"The selected test with the highest fault coverage is {test_name} with a complexity of {test_complexity}N.")
    print("The total test time is calculated using the formula:")
    print("Time = (RAM Size) * (Test Complexity) * (Time per Cycle)\n")
    print("Here is the equation with the specified values:")
    print(f"Time = {ram_size_bits} * {test_complexity} * {time_per_cycle_ns} ns")
    
    # Calculate intermediate result to show the total nanoseconds
    calculated_ns = ram_size_bits * test_complexity * time_per_cycle_ns
    print(f"Time = {int(calculated_ns)} ns")

    print("\nConverting the result to milliseconds:")
    print(f"Time = {int(calculated_ns)} ns / 1,000,000 ns per ms")
    print(f"Final Test Duration = {total_time_ms} ms")

# Execute the function to print the results.
calculate_ram_test_time()
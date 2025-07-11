def calculate_efficiency():
    """
    Calculates and prints the efficiency of the cell sorting experiment.
    """
    # The number of wells from which paired heavy and light chains were successfully recovered.
    successful_wells = 4

    # The total number of wells that were sorted with presumed cell doublets.
    total_sorted_wells = 95

    # Check to avoid division by zero, although not necessary with the given numbers.
    if total_sorted_wells == 0:
        print("Total sorted wells cannot be zero.")
        return

    # Calculate the efficiency percentage.
    efficiency = (successful_wells / total_sorted_wells) * 100

    # Print the equation with the numbers and the final result.
    print("Experimental Efficiency Calculation:")
    print(f"Equation: ({successful_wells} successful wells / {total_sorted_wells} total sorted wells) * 100")
    print(f"Result: {efficiency:.2f}%")
    print("\nThe success rate is extremely low, suggesting a fundamental challenge with the experimental assumptions.")

# Run the calculation.
calculate_efficiency()
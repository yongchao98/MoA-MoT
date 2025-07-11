def calculate_efficiency():
    """
    Calculates and explains the efficiency of the cell sorting experiment.
    """
    # Number of wells where paired heavy and light chains were successfully amplified
    successful_wells = 4

    # Total number of wells sorted with presumed doublets
    total_sorted_wells = 95

    # Calculate the experimental efficiency as a percentage
    efficiency = (successful_wells / total_sorted_wells) * 100

    print("--- Experimental Efficiency Calculation ---")
    print(f"Number of successful wells: {successful_wells}")
    print(f"Total number of sorted wells: {total_sorted_wells}")
    print("\nCalculation:")
    print(f"({successful_wells} / {total_sorted_wells}) * 100 = {efficiency:.2f}%")
    print(f"\nThe experimental efficiency was only {efficiency:.2f}%.")
    print("\n--- Analysis of Low Efficiency ---")
    print("The most fundamental reason for this low efficiency is the biological rarity of the event being studied.")

calculate_efficiency()

print("\n<<<B>>>")
def calculate_efficiency():
    """
    Calculates the efficiency of recovering paired antibody chains from sorted doublets.
    """
    # The number of wells where paired heavy and light chains were successfully amplified.
    successful_wells = 4

    # The total number of wells that were sorted with a putative doublet.
    total_sorted_wells = 95

    # The problem asks us to show the final equation with the numbers.
    print(f"Experimental Goal: Isolate paired antibody heavy and light chains from B cell:Tumor cell doublets.")
    print(f"Successful Outcomes: {successful_wells}")
    print(f"Total Attempts (Sorted Wells): {total_sorted_wells}")

    # Calculate the efficiency as a percentage.
    efficiency = (successful_wells / total_sorted_wells) * 100

    print("\nFinal Equation:")
    print(f"({successful_wells} / {total_sorted_wells}) * 100 = {efficiency:.2f}%")
    print(f"\nThe experimental efficiency is only {efficiency:.2f}%.")

calculate_efficiency()
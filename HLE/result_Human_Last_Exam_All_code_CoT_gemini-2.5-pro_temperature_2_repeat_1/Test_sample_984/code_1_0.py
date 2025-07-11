def calculate_efficiency(successful_events, total_events):
    """
    Calculates the percentage efficiency of an experimental outcome.

    Args:
        successful_events (int): The number of successful outcomes.
        total_events (int): The total number of attempts or trials.

    Returns:
        float: The calculated efficiency as a percentage.
    """
    if total_events == 0:
        return 0.0
    efficiency = (successful_events / total_events) * 100
    return efficiency

# --- Main part of the script ---
# Numbers from the problem description
successful_wells = 4
total_sorted_wells = 95

# Calculate the efficiency
efficiency = calculate_efficiency(successful_wells, total_sorted_wells)

# Print the analysis
print(f"Out of {total_sorted_wells} sorted wells, paired heavy and light chains were obtained from {successful_wells} wells.")
print(f"The equation for the experimental efficiency is: ({successful_wells} / {total_sorted_wells}) * 100%")
print(f"The calculated efficiency is: {efficiency:.2f}%")
print("\nThis very low efficiency highlights the rarity of the target event. The most plausible explanation is that true antigen-specific interactions are rare, and most sorted doublets were the result of random, non-specific cell association.")

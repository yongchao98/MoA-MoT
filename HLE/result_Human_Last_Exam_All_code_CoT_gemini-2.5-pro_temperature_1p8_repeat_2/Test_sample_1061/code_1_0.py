def simulate_isotope_change(start_year, end_year):
    """
    Simulates the effect of changing monsoon strength on the 13C isotope ratio in tree rings.

    This model makes the following assumptions:
    1. A "Monsoon Strength Index" ranges from 0 (very dry) to 1 (very wet).
    2. We model a hypothetical linear increase in monsoon strength from 0.4 to 0.75
       over the period 1886-1990.
    3. The 13C ratio decreases as monsoon strength increases, reflecting greater
       isotopic discrimination by the tree under less water stress.
    """
    print(f"Simulating the declining 13C ratio from {start_year} to {end_year} based on changing monsoon strength.\n")

    # Starting values
    base_13c_ratio = -25.0  # A typical baseline 13C ratio in permil (‰)
    start_monsoon_strength = 0.4
    end_monsoon_strength = 0.75

    # This factor determines how much the 13C ratio changes per unit of monsoon strength.
    # A negative value means the 13C ratio decreases as the monsoon gets stronger.
    sensitivity_factor = -3.0

    # --- Calculations for a few key years ---
    years_to_show = [start_year, start_year + (end_year - start_year) // 2, end_year]

    print(f"{'Year':<8} | {'Equation':<50} | {'Resulting 13C Ratio (‰)':<10}")
    print("-" * 80)

    for year in years_to_show:
        # Calculate the progress through the time period (from 0.0 to 1.0)
        progress = (year - start_year) / (end_year - start_year)

        # Calculate the monsoon strength for the given year
        current_monsoon_strength = start_monsoon_strength + (end_monsoon_strength - start_monsoon_strength) * progress

        # Calculate the change from the baseline 13C ratio
        change_in_13c = (current_monsoon_strength - start_monsoon_strength) * sensitivity_factor

        # Calculate the final 13C ratio
        final_13c_ratio = base_13c_ratio + change_in_13c

        # Print the full equation for clarity
        equation_str = (
            f"{base_13c_ratio:.2f} + ({current_monsoon_strength:.2f} - {start_monsoon_strength:.2f}) * {sensitivity_factor:.2f}"
        )
        print(f"{year:<8} | {equation_str:<50} | {final_13c_ratio:.3f}")


# Run the simulation for the period mentioned in the question
simulate_isotope_change(1886, 1990)

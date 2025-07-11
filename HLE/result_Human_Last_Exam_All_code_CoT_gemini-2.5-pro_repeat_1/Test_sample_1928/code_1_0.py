import math

def verify_statement_three():
    """
    This function verifies Statement III by comparing the cumulative duration of the
    three shortest U.S. presidential administrations to the orbital period of 16 Cygni Bb.
    """
    # Statement III: The cumulative duration of the three shortest U.S. presidential administrations
    # could fit within a local year at this location.

    # Orbital period of 16 Cygni Bb in Earth days (a "local year")
    # Source: NASA Exoplanet Archive
    orbital_period_16_cygni_bb = 799.5  # days

    # The three shortest U.S. presidential administrations and their durations in days
    shortest_administrations = {
        "William Henry Harrison": 31,
        "James A. Garfield": 199,
        "Zachary Taylor": 492
    }

    # Calculate the cumulative duration
    cumulative_duration = sum(shortest_administrations.values())

    print("Evaluating Statement III:")
    print(f"A local year on 16 Cygni Bb (its orbital period) is {orbital_period_16_cygni_bb} Earth days.")
    print("\nThe three shortest U.S. presidential administrations are:")
    for president, duration in shortest_administrations.items():
        print(f"- {president}: {duration} days")

    print(f"\nTheir cumulative duration is the sum of these periods:")
    calculation_str = " + ".join(map(str, shortest_administrations.values()))
    print(f"Equation: {calculation_str} = {cumulative_duration} days")

    print(f"\nComparing the two durations:")
    print(f"Cumulative administration duration ({cumulative_duration} days) vs. 16 Cygni Bb's year ({orbital_period_16_cygni_bb} days).")

    if cumulative_duration < orbital_period_16_cygni_bb:
        print("\nResult: The cumulative duration is SHORTER than a year on 16 Cygni Bb.")
        print("Therefore, Statement III is TRUE.")
    else:
        print("\nResult: The cumulative duration is LONGER than a year on 16 Cygni Bb.")
        print("Therefore, Statement III is FALSE.")

if __name__ == '__main__':
    verify_statement_three()
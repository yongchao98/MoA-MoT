import datetime

def solve_statement_iii():
    """
    This function verifies Statement III by comparing the orbital period of
    16 Cygni Bb with the combined duration of the three shortest U.S.
    presidential administrations.
    """
    # Orbital period of 16 Cygni Bb in Earth days.
    planet_orbital_period_days = 799.5

    # Durations of the three shortest U.S. presidential administrations in days.
    # 1. William Henry Harrison (Mar 4, 1841 – Apr 4, 1841)
    harrison_duration_days = 31
    # 2. James A. Garfield (Mar 4, 1881 – Sep 19, 1881)
    garfield_duration_days = 199
    # 3. Zachary Taylor (Mar 4, 1849 – Jul 9, 1850)
    taylor_duration_days = 492

    # Calculate the total duration of the three administrations.
    total_presidency_duration = harrison_duration_days + garfield_duration_days + taylor_duration_days

    # Compare the cumulative duration to the planet's orbital period.
    is_statement_true = total_presidency_duration < planet_orbital_period_days

    # Print the detailed analysis and result.
    print("--- Analysis of Statement III ---")
    print(f"A 'local year' on 16 Cygni Bb (its orbital period) is {planet_orbital_period_days} Earth days.")
    print("\nThe combined duration of the three shortest U.S. presidential administrations is calculated:")
    
    # The final equation with each number.
    print(f"Equation: {harrison_duration_days} days (Harrison) + {garfield_duration_days} days (Garfield) + {taylor_duration_days} days (Taylor) = {total_presidency_duration} days.")
    
    print(f"\nThis total of {total_presidency_duration} days is compared to the planet's year of {planet_orbital_period_days} days.")
    
    if is_statement_true:
        print("Result: The combined duration is SHORTER than a year on 16 Cygni Bb.")
        print("Conclusion: Statement III is TRUE.")
    else:
        print("Result: The combined duration is LONGER than a year on 16 Cygni Bb.")
        print("Conclusion: Statement III is FALSE.")

# Execute the function to print the analysis.
solve_statement_iii()
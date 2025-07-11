def check_statement_three():
    """
    This function verifies statement III by comparing the orbital period of
    16 Cygni Bb with the cumulative duration of the three shortest US
    presidential administrations.
    """
    # Orbital period of 16 Cygni Bb in Earth days
    orbital_period_days = 799.5

    # Durations of the three shortest U.S. presidential administrations in days
    harrison_term_days = 31
    garfield_term_days = 199
    taylor_term_days = 492

    # Calculate the cumulative duration
    cumulative_duration_days = harrison_term_days + garfield_term_days + taylor_term_days

    # Print the explanation and the numbers involved in the calculation
    print("Evaluating Statement III: Can the cumulative duration of the three shortest U.S. presidential administrations fit within a local year at 16 Cygni Bb?")
    print("-" * 50)
    print(f"A 'local year' (orbital period) on 16 Cygni Bb = {orbital_period_days} days.")
    print("\nDurations of the three shortest U.S. presidential administrations:")
    print(f"  - William Henry Harrison: {harrison_term_days} days")
    print(f"  - James A. Garfield: {garfield_term_days} days")
    print(f"  - Zachary Taylor: {taylor_term_days} days")
    
    # Print the equation
    print(f"\nSumming these terms: {harrison_term_days} + {garfield_term_days} + {taylor_term_days} = {cumulative_duration_days} days.")

    # Compare the two values and print the conclusion
    is_shorter = cumulative_duration_days < orbital_period_days
    print(f"\nComparison: Is {cumulative_duration_days} days < {orbital_period_days} days?")
    print(f"Result: {is_shorter}")
    
    if is_shorter:
        print("\nConclusion: The cumulative duration is shorter than the orbital period. Therefore, statement III is TRUE.")
    else:
        print("\nConclusion: The cumulative duration is not shorter than the orbital period. Therefore, statement III is FALSE.")

check_statement_three()
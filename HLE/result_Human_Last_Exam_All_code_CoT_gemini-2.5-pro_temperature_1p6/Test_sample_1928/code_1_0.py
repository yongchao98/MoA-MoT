def check_presidential_terms_vs_orbital_period():
    """
    Checks if the cumulative duration of the three shortest U.S. presidential
    administrations is less than the orbital period of 16 Cygni Bb.
    """
    # Orbital period of 16 Cygni Bb in Earth days
    orbital_period_16_cygni_bb = 799.5

    # Durations of the three shortest U.S. presidential administrations in days
    # 1. William Henry Harrison (March 4, 1841 - April 4, 1841)
    harrison_term = 31
    # 2. James A. Garfield (March 4, 1881 - September 19, 1881)
    garfield_term = 199
    # 3. Zachary Taylor (March 4, 1849 - July 9, 1850)
    taylor_term = 492

    # Calculate the sum of the durations
    sum_of_terms = harrison_term + garfield_term + taylor_term

    # Compare the sum to the orbital period
    is_shorter = sum_of_terms < orbital_period_16_cygni_bb

    # Print the results
    print("Evaluating Statement III:")
    print("-------------------------")
    print(f"Orbital Period of 16 Cygni Bb: {orbital_period_16_cygni_bb} days")
    print("\nThree Shortest U.S. Presidential Administrations:")
    print(f"1. William Henry Harrison: {harrison_term} days")
    print(f"2. James A. Garfield: {garfield_term} days")
    print(f"3. Zachary Taylor: {taylor_term} days")
    print("---")
    print(f"Equation: {harrison_term} + {garfield_term} + {taylor_term} = {sum_of_terms} days")
    print(f"Comparison: {sum_of_terms} days < {orbital_period_16_cygni_bb} days is {is_shorter}.")
    print("\nConclusion: The cumulative duration of the three shortest U.S. presidential administrations could indeed fit within a local year at 16 Cygni Bb. Statement III is true.")


check_presidential_terms_vs_orbital_period()

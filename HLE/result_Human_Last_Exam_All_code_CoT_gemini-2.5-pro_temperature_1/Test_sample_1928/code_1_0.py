import datetime

def verify_statement_three():
    """
    This function verifies Statement III by comparing the cumulative duration of the three
    shortest U.S. presidential administrations to the orbital period of 16 Cygni Bb.
    """
    # Statement III: The cumulative duration of the three shortest U.S. presidential administrations
    # could fit within a local year at this location.

    # A "local year" on 16 Cygni Bb is its orbital period.
    # Source: NASA Exoplanet Archive
    orbital_period_days = 803.88

    # The three shortest U.S. presidential administrations by term length.
    # 1. William Henry Harrison (died in office)
    harrison_term_days = 31
    # 2. James A. Garfield (assassinated)
    garfield_term_days = (datetime.date(1881, 9, 19) - datetime.date(1881, 3, 4)).days
    # 3. Zachary Taylor (died in office)
    taylor_term_days = (datetime.date(1850, 7, 9) - datetime.date(1849, 3, 4)).days

    # Calculate the cumulative duration of these administrations.
    cumulative_duration_days = harrison_term_days + garfield_term_days + taylor_term_days

    # Compare the cumulative duration to the planet's orbital period.
    is_statement_true = cumulative_duration_days < orbital_period_days

    # Print the detailed verification process.
    print("Verification of Statement III:")
    print(f"A local year on 16 Cygni Bb (its orbital period) is {orbital_period_days:.2f} Earth days.")
    print("\nThe three shortest U.S. presidential administrations are:")
    print(f"1. William Henry Harrison's term: {harrison_term_days} days")
    print(f"2. James A. Garfield's term: {garfield_term_days} days")
    print(f"3. Zachary Taylor's term: {taylor_term_days} days")

    print("\nThe final equation for the cumulative duration is:")
    print(f"{harrison_term_days} + {garfield_term_days} + {taylor_term_days} = {cumulative_duration_days} days")

    print("\nComparing the cumulative duration to the planet's local year:")
    print(f"Is {cumulative_duration_days} days < {orbital_period_days:.2f} days?")
    print(f"Result: {is_statement_true}")

    if is_statement_true:
        print("\nConclusion: The statement is TRUE as the cumulative duration is less than the orbital period.")
    else:
        print("\nConclusion: The statement is FALSE as the cumulative duration is not less than the orbital period.")

verify_statement_three()
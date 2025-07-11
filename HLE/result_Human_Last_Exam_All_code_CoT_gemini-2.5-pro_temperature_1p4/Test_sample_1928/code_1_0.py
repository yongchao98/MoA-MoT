#
# This script verifies Statement III by comparing the orbital period of 16 Cygni Bb
# with the cumulative duration of the three shortest U.S. presidential administrations.
#

def verify_statement_three():
    """
    Performs the calculation to verify Statement III and prints the results.
    """
    # Orbital period of 16 Cygni Bb in Earth days (a "local year")
    period_16_cygni_bb = 799.5

    # Durations of the three shortest U.S. presidential administrations in days
    harrison_days = 31
    garfield_days = 199
    taylor_days = 492

    # Calculate the cumulative duration
    total_short_presidencies = harrison_days + garfield_days + taylor_days

    # Print the findings
    print("Verifying Statement III:")
    print(f"A local year on 16 Cygni Bb (its orbital period) is {period_16_cygni_bb} days.")
    print("\nThe three shortest U.S. presidential administrations are:")
    print(f"- William Henry Harrison: {harrison_days} days")
    print(f"- James A. Garfield: {garfield_days} days")
    print(f"- Zachary Taylor: {taylor_days} days")
    print("\nCalculating the cumulative duration:")
    # The prompt requests that the final equation with all numbers be printed
    print(f"{harrison_days} + {garfield_days} + {taylor_days} = {total_short_presidencies} days.")

    # Compare the two durations and print the conclusion
    if total_short_presidencies < period_16_cygni_bb:
        print(f"\nConclusion: Since {total_short_presidencies} days is less than {period_16_cygni_bb} days, the statement is TRUE.")
    else:
        print(f"\nConclusion: Since {total_short_presidencies} days is not less than {period_16_cygni_bb} days, the statement is FALSE.")

if __name__ == "__main__":
    verify_statement_three()

import sys

def find_rarest_noble_gas_percentage():
    """
    This script identifies the rarest noble gas and calculates its approximate
    percentage in terrestrial matter.
    """
    # Step 1: Identify the rarest noble gas and its abundance.
    # The rarest naturally occurring noble gas on Earth is Radon (Rn).
    # This is due to its high radioactivity and short half-life.
    # Its approximate mass fraction in the atmosphere is 1 part per 10^21.
    gas_name = "Radon"
    numerator = 1.0
    denominator = 1e21  # 1e21 is scientific notation for 1 followed by 21 zeros.

    # Step 2: Calculate the percentage.
    # The formula is (fraction) * 100.
    percentage = (numerator / denominator) * 100

    # Step 3: Print the results, showing the numbers in the final equation.
    print(f"The rarest noble gas on Earth (in 2002 and today) is {gas_name}.")
    print(f"Its approximate mass fraction is {int(numerator)} part in {denominator:.0e}.")
    print("\nThe calculation to find the percentage is:")
    # sys.stdout.write is used to prevent the f-string from evaluating the equation early.
    sys.stdout.write(f"({int(numerator)} / {denominator:.0e}) * 100 = {percentage:.1e}%\n")
    print(f"\nSo, the percentage is approximately {percentage:.1e}%.")

if __name__ == "__main__":
    find_rarest_noble_gas_percentage()
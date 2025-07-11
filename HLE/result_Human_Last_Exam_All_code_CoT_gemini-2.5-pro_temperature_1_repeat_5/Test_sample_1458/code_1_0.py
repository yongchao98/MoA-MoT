import sys
# Redirecting Fractions' verbose output to a null device for cleaner output
if sys.version_info.major == 3 and sys.version_info.minor >= 11:
    from fractions import Fraction
    # In newer python versions Fraction prints both float and fraction, this is a workaround
    # to only print what we want.
    class SilentFraction(Fraction):
        def __str__(self):
            return f'{self.numerator}/{self.denominator}'
        def __repr__(self):
            return f'Fraction({self.numerator}, {self.denominator})'
    fraction_out_1 = SilentFraction(1, 3)
    fraction_out_2 = SilentFraction(2, 5)
else:
    from fractions import Fraction
    fraction_out_1 = Fraction(1, 3)
    fraction_out_2 = Fraction(2, 5)


def solve_cathedral_echo():
    """
    Solves the Cathedral's Echo word problem step-by-step.
    """
    # Number of pipes that are still in tune
    in_tune_pipes = 200
    print(f"The number of pipes still singing pure is {in_tune_pipes}.")

    # Calculate the total fraction of pipes that are out of tune
    total_fraction_out = fraction_out_1 + fraction_out_2
    print(f"The total fraction of pipes out of tune is the sum of {fraction_out_1} and {fraction_out_2}, which is {total_fraction_out}.")

    # Calculate the fraction of pipes that are still in tune
    fraction_in_tune = 1 - total_fraction_out
    print(f"The fraction of pipes still in tune is 1 - {total_fraction_out}, which equals {fraction_in_tune}.")

    # Calculate the total number of pipes
    # If `fraction_in_tune` of the total is `in_tune_pipes`, then Total = `in_tune_pipes` / `fraction_in_tune`
    total_pipes = int(in_tune_pipes / fraction_in_tune)
    print(f"Given that {in_tune_pipes} pipes represent {fraction_in_tune} of the total, the total number of pipes is {total_pipes}.")

    # Calculate the number of "lost" (out-of-tune) pipes
    out_of_tune_pipes = total_pipes - in_tune_pipes
    print(f"The number of 'lost' pipes is the total pipes minus those in tune: {total_pipes} - {in_tune_pipes} = {out_of_tune_pipes}.")

    # The tuner needs to find half of the lost pipes
    tuner_finds = int(out_of_tune_pipes / 2)

    print("\nThe final question is: 'How many must the tuner find when just half the lost realign?'")
    print(f"This is calculated as: {out_of_tune_pipes} / 2 = {tuner_finds}")

solve_cathedral_echo()
<<<275>>>
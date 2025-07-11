import fractions

def calculate_and_print_difference():
    """
    Calculates and prints the solution to the Yahtzee probability problem.
    """
    # --- Strategy A: Aim for five '1's ---

    # Probability of rolling a '1' is p(1) = 2^(-1) = 1/2
    p1 = fractions.Fraction(1, 2)

    # Probability of a single re-rolled die becoming a '1' in at most two tries
    q1 = p1 + (1 - p1) * p1

    # Probability of turning four dice into '1's
    prob_A = q1 ** 4

    # --- Strategy B: Aim for five '3's ---

    # Probability of rolling a '3' is p(3) = 2^(-3) = 1/8
    p3 = fractions.Fraction(1, 8)

    # Probability of a single re-rolled die becoming a '3' in at most two tries
    q3 = p3 + (1 - p3) * p3

    # Probability of turning two dice into '3's
    prob_B = q3 ** 2

    # --- Calculate and Print the Difference ---

    # Find the difference between the two probabilities
    diff = prob_A - prob_B
    prob_A_common_denom = prob_A.limit_denominator(diff.denominator)

    print("--- Strategy A: Keep the '1' and aim for five '1's ---")
    print(f"The probability of a die becoming a '1' in two rolls is ({p1.numerator}/{p1.denominator}) + (1 - {p1.numerator}/{p1.denominator}) * ({p1.numerator}/{p1.denominator}) = {q1.numerator}/{q1.denominator}.")
    print(f"The total probability for this strategy is P(A) = ({q1.numerator}/{q1.denominator})^4 = {prob_A.numerator}/{prob_A.denominator}.\n")

    print("--- Strategy B: Keep the three '3's and aim for five '3's ---")
    print(f"The probability of a die becoming a '3' in two rolls is ({p3.numerator}/{p3.denominator}) + (1 - {p3.numerator}/{p3.denominator}) * ({p3.numerator}/{p3.denominator}) = {q3.numerator}/{q3.denominator}.")
    print(f"The total probability for this strategy is P(B) = ({q3.numerator}/{q3.denominator})^2 = {prob_B.numerator}/{prob_B.denominator}.\n")

    print("--- Final Difference Calculation ---")
    print("The difference in probabilities between the two strategies is P(A) - P(B).")
    print(f"Difference = {prob_A.numerator}/{prob_A.denominator} - {prob_B.numerator}/{prob_B.denominator}")
    print(f"Difference = {prob_A_common_denom.numerator}/{prob_A_common_denom.denominator} - {prob_B.numerator}/{prob_B.denominator} = {diff.numerator}/{diff.denominator}")

if __name__ == '__main__':
    calculate_and_print_difference()
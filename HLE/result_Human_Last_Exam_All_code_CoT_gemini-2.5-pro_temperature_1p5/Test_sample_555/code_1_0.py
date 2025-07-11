from fractions import Fraction

def find_minimal_expected_rolls():
    """
    This function calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit from 0 to 9.

    It analyzes strategies based on grouping rolls into blocks of size k,
    and finds the k that minimizes the expected number of rolls.
    """
    min_expected_rolls = float('inf')
    best_k = -1
    min_fraction = None

    print("Analyzing strategies based on rolling the die k times in each trial.")
    print("-" * 60)

    # We test for k from 2 upwards. We can stop when k itself is greater
    # than the current minimum expectation, since E_k > k.
    for k in range(2, 10):
        power_of_7 = 7**k
        
        # This is the largest multiple of 10 <= 7^k
        num_usable_outcomes = 10 * (power_of_7 // 10)

        if num_usable_outcomes == 0:
            # This case happens for k=1, but the loop starts at k=2.
            # It's included for completeness.
            continue

        # Calculate the expected number of rolls E_k as a fraction
        # E_k = k / p_k = k / (num_usable_outcomes / power_of_7)
        expected_rolls_k = Fraction(k * power_of_7, num_usable_outcomes)

        print(f"For a strategy with k = {k} rolls per trial:")
        print(f"  - Total outcomes per trial: 7^{k} = {power_of_7}")
        print(f"  - Usable outcomes for 10 digits: 10 * floor({power_of_7}/10) = {num_usable_outcomes}")
        print(f"  - Expected rolls E_{k} = ({k} * {power_of_7}) / {num_usable_outcomes} = {expected_rolls_k}")
        print(f"  - Decimal value: E_{k} ≈ {float(expected_rolls_k):.4f}")
        print("-" * 25)

        if expected_rolls_k < min_expected_rolls:
            min_expected_rolls = expected_rolls_k
            best_k = k
        
        # Optimization: If k becomes larger than the best E_k found so far,
        # no future k' can be better, since E_k' > k'.
        if k > min_expected_rolls:
            break

    print("\nConclusion:")
    print(f"The minimal expected value is achieved with a strategy using k = {best_k} rolls per trial.")
    
    # Re-calculate and display the equation for the best k
    k = best_k
    power_of_7 = 7**k
    denominator = 10 * (power_of_7 // 10)
    numerator = k * power_of_7
    final_fraction = Fraction(numerator, denominator)
    
    print(f"\nThe minimal expected number of rolls is calculated as:")
    print(f"E_{k} = ({k} * 7^{k}) / (10 * floor(7^{k}/10))")
    print(f"E_{k} = ({k} * {power_of_7}) / (10 * {power_of_7//10})")
    print(f"E_{k} = {numerator} / {denominator} = {final_fraction.numerator}/{final_fraction.denominator}")

if __name__ == '__main__':
    find_minimal_expected_rolls()
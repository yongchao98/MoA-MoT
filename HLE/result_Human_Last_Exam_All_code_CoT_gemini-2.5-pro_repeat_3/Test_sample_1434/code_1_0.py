import math

def calculate_series_sum():
    """
    This function calculates the sum based on the analysis of the problem.
    The sum is split into two parts: one for odd n and one for even n.
    Both parts form a geometric series with ratio r = 1/8.
    """

    # The sum for odd n is Sum_odd = 4 * sum_{j=1 to inf} (1/8)^j
    # The sum of the geometric series is (1/8) / (1 - 1/8) = 1/7
    sum_odd_numerator = 4
    sum_odd_denominator = 7
    
    # The sum for even n is Sum_even = 2 * sum_{j=1 to inf} (1/8)^j
    sum_even_numerator = 2
    sum_even_denominator = 7

    # Total sum is Sum_odd + Sum_even
    total_sum_numerator = sum_odd_numerator + sum_even_numerator
    total_sum_denominator = sum_odd_denominator

    print("The calculation is split into a sum over odd n and a sum over even n.")
    print("Sum over odd n:")
    print(f"Numerator: {sum_odd_numerator}")
    print(f"Denominator: {sum_odd_denominator}")

    print("\nSum over even n:")
    print(f"Numerator: {sum_even_numerator}")
    print(f"Denominator: {sum_even_denominator}")
    
    print("\nThe final sum is the addition of these two parts.")
    print("Final equation:")
    # As requested, printing each number in the final equation
    print(f"{sum_odd_numerator}/{sum_odd_denominator} + {sum_even_numerator}/{sum_even_denominator} = {total_sum_numerator}/{total_sum_denominator}")

calculate_series_sum()

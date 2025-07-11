import math

def compute_series_sum():
    """
    This function computes the value of the series S = sum_{i=1 to inf} 1 / (2^i * (2^i - 1)).
    The series converges very quickly, so summing a finite number of terms gives a good approximation.
    """
    series_sum = 0.0
    # Summing up to 60 terms is more than enough for standard double precision.
    num_terms = 60
    for i in range(1, num_terms + 1):
        power_of_2 = 1 << i  # Efficient way to calculate 2**i for integers
        term = 1.0 / (power_of_2 * (power_of_2 - 1))
        series_sum += term
    return series_sum

# The total sum is 2 * T - 1/3, where T is the series sum computed above.
series_t_value = compute_series_sum()
factor = 2
subtrahend = 1/3
result = factor * series_t_value - subtrahend

# As requested, here we output the numbers in the final equation.
print("The final sum is calculated based on the equation: 2 * T - 1/3")
print(f"The value of T, which is sum(1/(2^i(2^i-1))), is approximated as: {series_t_value}")
print(f"Final calculation: {factor} * {series_t_value} - {subtrahend} = {result}")

# Per instructions, the final answer in the specified format is below this block.
# print(f"<<<{result}>>>")
import math

def calculate_av_n_k(n):
    """
    Calculates the number of 1324-avoiding permutations of length n with 3 inversions.
    The formula used is (n^3 - 6*n^2 + 11*n - 4) / 2.
    """
    if n < 3:
        return 0

    # Calculate each term of the numerator
    term1 = n**3
    term2 = -6 * n**2
    term3 = 11 * n
    term4 = -4

    numerator = term1 + term2 + term3 + term4
    result = numerator // 2

    # Output the equation with the calculated numbers
    print(f"({term1} + ({term2}) + {term3} + ({term4})) / 2 = {result}")

    return result

# The value of n for the problem
n = 333
print(f"Calculating av_{n}^3(1324):")

# Call the function to compute and print the result
final_answer = calculate_av_n_k(n)
print(f"The value of av_{333}^3(1324) is {final_answer}")

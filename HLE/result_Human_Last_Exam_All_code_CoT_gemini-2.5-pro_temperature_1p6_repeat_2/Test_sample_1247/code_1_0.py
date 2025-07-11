import math

def solve():
    """
    Calculates av_n^k(1324) for n=333 and k=3 using a known formula.
    av_n^3(1324) = (n^3 + 6n^2 - 43n + 54) / 6 for n >= 5.
    """
    n = 333

    # Calculate each term of the polynomial
    term1 = n**3
    term2 = 6 * (n**2)
    term3 = 43 * n
    term4 = 54

    # Calculate the numerator
    numerator = term1 + term2 - term3 + term4

    # The result is the numerator divided by 6
    result = numerator // 6

    # Print the equation with the evaluated numbers and the final answer
    print(f"The number of 1324-avoiding permutations of length {n} with 3 inversions is calculated using the formula:")
    print(f"av_{n}^3(1324) = ({n}^3 + 6*{n}^2 - 43*{n} + 54) / 6")
    print(f"av_{n}^3(1324) = ({term1} + {term2} - {term3} + {term4}) / 6")
    print(f"av_{n}^3(1324) = {numerator} / 6")
    print(f"av_{n}^3(1324) = {result}")

solve()

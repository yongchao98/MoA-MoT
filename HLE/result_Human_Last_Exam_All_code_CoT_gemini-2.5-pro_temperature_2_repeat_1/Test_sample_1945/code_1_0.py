import math

def get_limit():
    """
    Calculates the theoretical limit for Pr(n) as n -> infinity.
    The limit is given by the equation: L = 3 / pi^2
    """
    numerator = 3
    denominator = math.pi**2
    limit_value = numerator / denominator

    print("The necessary and sufficient condition for (p,q) to be a good pair is: p + q <= n + 1")
    print("\nThe limit of Pr(n) as n -> infinity is given by the equation: L = 3 / pi^2")
    print("The numbers in the final equation are:")
    print(f"Numerator: {numerator}")
    print(f"Denominator (pi^2): {denominator}")
    print(f"The exact value of the limit is: {limit_value}")

if __name__ == '__main__':
    get_limit()
import math

def solve():
    """
    Calculates the integer part of 10^4 * L, where L is the limit derived from the problem.
    L = 2/ln(alpha) + 2/ln(beta)
    """

    # The formula for the limit L
    # L = 2 / ln((13 + sqrt(165)) / 2) + 2 / ln((5 + sqrt(21)) / 2)

    # Outputting the numbers from the equation as requested
    print("The final calculation is based on the formula: 10^4 * (2 / ln((13 + sqrt(165)) / 2) + 2 / ln((5 + sqrt(21)) / 2))")
    print("The numbers in the equation are: 10000, 2, 13, 165, 2, 2, 5, 21, 2")

    # Define the growth rates alpha and beta
    alpha = (13 + math.sqrt(165)) / 2
    beta = (5 + math.sqrt(21)) / 2

    # Calculate the natural logarithms
    log_alpha = math.log(alpha)
    log_beta = math.log(beta)

    # Calculate the limit L
    L = 2 / log_alpha + 2 / log_beta

    # Calculate the final value
    result = 10**4 * L

    # Print the integer part of the result
    print(f"\nThe value of alpha is: {alpha}")
    print(f"The value of beta is: {beta}")
    print(f"The value of ln(alpha) is: {log_alpha}")
    print(f"The value of ln(beta) is: {log_beta}")
    print(f"The limit L is: {L}")
    print(f"The value of 10^4 * L is: {result}")
    print(f"\nThe integer part of the final answer is:")
    print(int(result))

solve()
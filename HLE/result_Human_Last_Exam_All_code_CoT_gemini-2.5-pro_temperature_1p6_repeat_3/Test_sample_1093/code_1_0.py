import math

def calculate_gauss_game_average():
    """
    Calculates the value that the average of the recorded numbers in Gauss's game
    converges to. This value is Euler's number, e.
    The function demonstrates this by computing e from its series expansion:
    e = 1/0! + 1/1! + 1/2! + 1/3! + ...
    """
    # The number of terms to sum for a high-precision approximation of e.
    # 18 terms is sufficient for standard float precision.
    num_terms = 18
    e_approximation = 0.0

    print("The average number of throws converges to Euler's number, e.")
    print("This can be calculated using the series expansion: e = 1/0! + 1/1! + 1/2! + 1/3! + ...")
    print("\nCalculating the value by summing each number in this equation:\n")

    for n in range(num_terms):
        term_value = 1.0 / math.factorial(n)
        e_approximation += term_value
        print(f"Term (1/{n}!): {term_value}")

    print("\n-------------------------------------------")
    print(f"The final average the numbers converge to is the sum of the series: {e_approximation}")
    print(f"(This is Euler's number, e, which is approximately 2.71828)")

calculate_gauss_game_average()
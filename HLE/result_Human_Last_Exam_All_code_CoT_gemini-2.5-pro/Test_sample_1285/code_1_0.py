def calculate_max_roots():
    """
    This function calculates the maximum number of roots for t=5
    based on the formula t*(t-1)/2 and prints the steps.
    """
    # Set the value of t for part (b)
    t = 5

    # The formula for the maximum number of roots is t * (t - 1) / 2
    
    # Calculate the numerator of the formula
    numerator = t * (t - 1)
    
    # The denominator of the formula
    denominator = 2
    
    # Calculate the final result
    max_roots = numerator // denominator

    # Print the equation for t=5, showing each number
    print(f"The calculation for t = {t} is as follows:")
    print(f"Maximum number of roots = ({t} * ({t} - 1)) / {denominator}")
    print(f"                          = {numerator} / {denominator}")
    print(f"                          = {max_roots}")

# Execute the function
calculate_max_roots()
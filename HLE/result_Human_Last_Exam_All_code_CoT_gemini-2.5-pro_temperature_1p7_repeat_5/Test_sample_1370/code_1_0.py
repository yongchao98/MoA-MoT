import math

def solve_dice_independence():
    """
    Calculates the largest possible number of mutually independent events
    for rolling 100 regular 6-sided dice.
    """
    num_dice = 100
    num_sides = 6

    # The size of the sample space N is num_sides raised to the power of num_dice.
    # According to a theorem in probability, the maximum number of mutually
    # independent events in an equiprobable sample space of size N is the
    # sum of the exponents in the prime factorization of N.

    # Prime factorization of N = 6^100 = (2 * 3)^100 = 2^100 * 3^100
    
    # The prime factors of the base (6) are 2 and 3.
    # The exponents in the factorization of N are therefore 100 for the prime 2,
    # and 100 for the prime 3.
    
    exponent_of_2 = 1 * num_dice
    exponent_of_3 = 1 * num_dice
    
    print(f"The size of the sample space is N = {num_sides}^{num_dice}.")
    print(f"The prime factorization of N is 2^{exponent_of_2} * 3^{exponent_of_3}.")
    
    # The maximum number of mutually independent events is the sum of the exponents.
    max_m = exponent_of_2 + exponent_of_3
    
    print("\nThe largest possible value of m is the sum of these exponents.")
    print("The final calculation is:")
    
    # Printing the numbers in the final equation as requested.
    print(f"{exponent_of_2} + {exponent_of_3} = {max_m}")
    
solve_dice_independence()
<<<200>>>
import math

def solve():
    """
    This function determines the first possible year the equation could have been written,
    based on the history of the symbols used.
    """
    
    # Year of introduction for various symbols:
    # Equals sign (=): 1557
    # Plus-minus sign (±): 1628
    # Newton's/Leibniz's calculus notation (ṡ, s̈, d/ds): c. 1684
    # At sign (@) in a technical/digital context: 1963 (Inclusion in ASCII)
    
    # The latest date is the limiting factor. The use of '@' as a variable is
    # a modern phenomenon, tied to its inclusion in standard character sets for computers.
    first_possible_year = 1963
    
    # Round the year to the nearest 10.
    # We can achieve this by dividing by 10, rounding, and then multiplying by 10.
    rounded_year = round(first_possible_year / 10) * 10
    
    print(f"The analysis of symbols points to the introduction of the '@' symbol in the ASCII standard as the limiting factor.")
    print(f"First possible year of publication: {first_possible_year}")
    print(f"Year rounded to the nearest 10 years: {rounded_year}")
    
    # The numbers present in the original equation are 1 and 2.
    print(f"The numbers in the equation are: 1, 2")

solve()
<<<1960>>>
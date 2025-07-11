import math

def solve():
    """
    Determines the earliest possible year an equation could have been written
    based on the history of its symbols.
    """
    
    # The prompt requests to output the numbers present in the equation.
    print("The numbers in the equation are 1 and 2.")
    
    # Step 1: Identify the symbol that was most recently adopted for technical use.
    # The equation uses '=', '±', and derivative notations from the 17th century.
    # The symbol '@' used as a variable is the most modern element.
    # Its use as a general technical symbol is tied to its inclusion in computer character sets.
    symbol_name = "'@' (at sign)"
    
    # Step 2: Determine the earliest year of its adoption in a technical context.
    # The ASCII character set, which included '@', was first published in 1963.
    # This is the earliest plausible year for an author to use it as a variable in a paper.
    earliest_year = 1963
    
    print(f"The limiting symbol is the {symbol_name} used as a mathematical variable.")
    print(f"Its first appearance in a standard technical character set (ASCII) was in {earliest_year}.")
    
    # Step 3: Round the year to the nearest 10 years.
    # We can round a number n to the nearest 10 by calculating round(n / 10) * 10.
    rounded_year = round(earliest_year / 10) * 10
    
    print(f"The earliest possible year is {earliest_year}.")
    print(f"Rounding {earliest_year} to the nearest 10 years gives {rounded_year}.")

solve()

print("<<<1960>>>")
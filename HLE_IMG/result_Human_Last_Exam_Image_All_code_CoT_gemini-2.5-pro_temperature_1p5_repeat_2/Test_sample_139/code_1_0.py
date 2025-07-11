import math

def solve():
    """
    This function determines the first possible year an equation could have been written,
    based on the historical introduction of its symbols, and rounds it to the nearest decade.
    """
    # Research shows the symbols in the equation and their first appearance dates are:
    # - Minus sign (-): 1489
    # - Equals sign (=): 1557
    # - Plus-minus sign (±): 1631
    # - Newton's dot notation (ṡ) for derivatives: ~1665
    # - Leibniz's notation (d/ds) for derivatives: ~1675
    # - At sign (@) used as a variable in a technical/mathematical context.
    #   While the symbol is old, its use in this manner is modern. A conservative
    #   estimate for its first possible use in a published paper is after it became
    #   widely known in the technical community via email.
    
    # The year the @ symbol was adopted for email, making it a familiar technical symbol.
    critical_year = 1971
    
    # The problem asks to round the year to the nearest 10 years.
    # In Python, round(number, -1) rounds to the nearest 10.
    # For example, round(1971, -1) becomes 1970.
    rounded_year = int(round(critical_year, -1))
    
    print(f"The limiting symbol is the '@' sign, first used in a widespread technical context around {critical_year}.")
    print(f"Rounding {critical_year} to the nearest 10 years gives {rounded_year}.")

solve()
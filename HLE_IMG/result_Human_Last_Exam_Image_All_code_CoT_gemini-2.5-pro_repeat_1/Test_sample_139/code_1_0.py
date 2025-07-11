import math

def solve_and_print():
    """
    This function determines the earliest possible year for the equation based on its symbols,
    and then rounds that year to the nearest 10.
    """
    # The numbers present in the equation and the accompanying text are 1, 2, and 0.
    # From the equation: @ = r / (1 - s^2 ± 2s) -> Numbers are 1 and 2.
    # From the text: d@/ds = 0 -> Number is 0.
    num1 = 1
    num2 = 2
    num0 = 0
    print(f"The numbers found in the equation are: {num1}, {num2}, {num0}")

    # The limiting symbol is '@'. Its widespread technical use began with email on ARPANET.
    # This event, initiated by Ray Tomlinson, occurred in 1971.
    # This is the most plausible "first possible year" for this symbol to appear in a published paper.
    first_possible_year = 1971

    # Round the year to the nearest 10.
    # We can use the round() function with a negative precision for this.
    rounded_year = round(first_possible_year, -1)

    print(f"The first possible year based on the symbols is {first_possible_year}.")
    print(f"Rounding to the nearest 10 years, the final answer is: {rounded_year}")

solve_and_print()
<<<1970>>>
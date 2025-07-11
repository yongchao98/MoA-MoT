def solve_haiku_riddle():
    """
    This function solves the riddle presented in the haiku.
    1. Extracts numbers from the haiku's text.
    2. Forms an equation based on the haiku's wording.
    3. Calculates the result and prints the full equation.
    The final answer is derived by interpreting this result.
    """
    # "An August tempest" -> August is the 8th month.
    num_august = 8

    # "Twice fifteen" -> 2 * 15
    num_twice_fifteen = 30

    # "A divine one" -> 1
    num_one = 1

    # The words suggest the equation:
    # "tempest" (destruction/loss) and "yields" (gives up) imply subtraction
    # from the main number, 30.
    # The equation is 30 - 8 - 1.
    result = num_twice_fifteen - num_august - num_one

    # The prompt requires printing each number in the final equation.
    # The result, 21, corresponds to the 21st letter of the alphabet: U.
    # The title "The Bays" is the question, and the answer is "U"
    # because bays are often U-shaped.
    print(f"The equation from the haiku is: {num_twice_fifteen} - {num_august} - {num_one} = {result}")

solve_haiku_riddle()
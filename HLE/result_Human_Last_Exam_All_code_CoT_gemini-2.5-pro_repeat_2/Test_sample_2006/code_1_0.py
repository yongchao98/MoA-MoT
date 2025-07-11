def solve_haiku_riddle():
    """
    This function solves the riddle presented in the haiku.
    It extracts numbers and operations from the poem, calculates the result,
    and converts it to the corresponding letter of the alphabet.
    """

    # Numbers derived from the haiku
    # "Twice fifteen brings winds of change" -> 2 * 15 = 30
    val_from_line2 = 30
    
    # "An August tempest" -> August is the 8th month. "Tempest" implies subtraction.
    val_from_line1 = 8
    
    # "A divine one yields" -> The number is 1. "Yields" implies subtraction.
    val_from_line3 = 1

    # The haiku describes the equation: 30 - 8 - 1
    result = val_from_line2 - val_from_line1 - val_from_line3

    # The question "in alphabetical order what is the answer?" suggests
    # converting the numerical result to a letter (1=A, 2=B, etc.).
    # We can get the character by adding the result to the ASCII value of 'A' - 1.
    # ASCII 'A' is 65. So for result 1, 1+64='A'. For 2, 2+64='B'.
    final_answer_char = chr(result + 64)

    print("The riddle in the haiku describes a mathematical equation.")
    print("Here is the breakdown:")
    print(f"'Twice fifteen' provides the number: {val_from_line2}")
    print(f"'August' provides the number: {val_from_line1}")
    print(f"'one' provides the number: {val_from_line3}")
    print("\nThe words 'tempest' and 'yields' imply subtraction.")
    print("\nThe resulting equation is:")
    print(f"{val_from_line2} - {val_from_line1} - {val_from_line3} = {result}")

    print(f"\nThe {result}st letter in alphabetical order is the answer.")
    print(f"The final answer is: {final_answer_char}")

solve_haiku_riddle()
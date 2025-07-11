import string

def solve_haiku_riddle():
    """
    This function solves the riddle presented in the haiku.
    It deciphers the numbers and operations from the text,
    calculates the result, and finds the corresponding letter of the alphabet.
    """
    
    # 1. Extract numbers from the haiku's phrases.
    # "Twice fifteen brings winds of change"
    num1 = 2 * 15
    # "An August tempest"
    num2 = 8
    # "A divine one yields"
    num3 = 1

    # 2. Perform the calculation. The wording "brings winds of change"
    # suggests subtraction from the primary number, 30.
    result = num1 - num2 - num3

    # 3. Find the corresponding letter in the alphabet.
    alphabet = string.ascii_uppercase
    # The result is a 1-based index for the alphabet.
    final_answer = alphabet[result - 1]

    # 4. Print the reasoning and the equation.
    print("The riddle can be solved by turning the phrases into a mathematical equation.")
    print("The numbers derived from the haiku are:")
    print(f"- 'Twice fifteen': {num1}")
    print(f"- 'An August tempest': {num2}")
    print(f"- 'A divine one': {num3}")
    print("\nThe equation formed is:")
    print(f"{num1} - {num2} - {num3} = {result}")
    
    print(f"\nThe {result}st letter of the alphabet is the answer: {final_answer}")
    print("The answer 'U' relates to the title 'The Bays', as bays are often U-shaped.")

solve_haiku_riddle()
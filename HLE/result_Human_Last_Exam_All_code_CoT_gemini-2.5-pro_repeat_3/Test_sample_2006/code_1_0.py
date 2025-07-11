import string

def solve_riddle():
    """
    Solves the haiku riddle by extracting numbers, performing a calculation,
    and finding the corresponding letter in the alphabet.
    """
    # Step 1: Extract numbers from the haiku's phrases
    # "Twice fifteen brings winds of change" -> 2 * 15 = 30
    twice_fifteen = 30
    
    # "An August tempest" -> August is the 8th month
    august_tempest = 8
    
    # "A divine one yields" -> 1
    divine_one = 1
    
    # Step 2 & 3: Perform the calculation
    # The central line (30) is the minuend, and the other two numbers are subtracted.
    result_number = twice_fifteen - august_tempest - divine_one
    
    # Step 4: Convert the number to its corresponding letter
    # A=1, B=2, ..., Z=26. The 21st letter is 'U'.
    alphabet = string.ascii_uppercase
    final_answer = alphabet[result_number - 1]
    
    # Print the equation and the final result
    print("The equation derived from the haiku is:")
    print(f"{twice_fifteen} - {august_tempest} - {divine_one} = {result_number}")
    print(f"The {result_number}st letter of the alphabet is the answer.")
    print(f"The answer to 'The Bays?' is: {final_answer}")

solve_riddle()
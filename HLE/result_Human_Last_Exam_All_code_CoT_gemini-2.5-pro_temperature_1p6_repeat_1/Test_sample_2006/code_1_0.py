import string

def solve_riddle():
    """
    This function solves the haiku riddle step-by-step.
    """
    print("Step 1: Decoding each line of the haiku into a number.")

    # "An August tempest" -> August is the 8th month.
    num1 = 8
    print(f"   - 'An August tempest': August is the 8th month, which gives the number {num1}.")

    # "Twice fifteen brings winds of change" -> 2 * 15
    num2 = 30
    print(f"   - 'Twice fifteen...': This translates to 2 * 15, which gives the number {num2}.")

    # "A divine one yields" -> The number one.
    num3 = 1
    print(f"   - 'A divine one yields': This explicitly gives the number {num3}.")

    print("\nStep 2: Summing the numbers from the haiku.")
    total = num1 + num2 + num3
    # The final equation is printed with each number.
    print(f"The final equation is: {num1} + {num2} + {num3} = {total}")

    print("\nStep 3: Finding the answer 'in alphabetical order'.")
    print("This means the answer is the letter of the alphabet corresponding to the sum (39).")
    print("Since there are only 26 letters, we wrap around using modular arithmetic.")

    # To find the Nth letter of the alphabet (A=1, B=2...), we can use the formula:
    # index = (N - 1) % 26
    # For N=39, index = (39 - 1) % 26 = 38 % 26 = 12.
    # The letter at index 12 in the alphabet (starting from 0) is 'M'.
    alphabet = string.ascii_uppercase
    final_letter_index = (total - 1) % 26
    final_answer = alphabet[final_letter_index]

    print(f"The 39th letter is the same as the {final_letter_index + 1}th letter, which is:")
    print(f"The final answer is: {final_answer}")

solve_riddle()
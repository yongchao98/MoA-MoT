import sys
import io

# A buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_riddle():
    """
    Solves the riddle by revealing the word from Joseph Heller's 'Catch-22'.
    """
    answer_word = "illiterate"
    
    print("The riddle is a reference to a passage in Joseph Heller's 1961 novel, 'Catch-22'.")
    print("In the book, a character tells a fabricated story about Pope Paul II to make a point.")
    print("The final word of the punchline, 'X', is what the riddle is asking for.")
    print("\nHere is the decoded word, with each character's ASCII value shown as an equation:\n")

    equations = {
        'i': "100 + 5",
        'l': "100 + 8",
        't': "110 + 6",
        'e': "100 + 1",
        'r': "110 + 4",
        'a': "90 + 7"
    }
    
    final_equation_parts = []
    for char in answer_word:
        equation_str = equations[char]
        numbers = [int(n) for n in equation_str.replace('+', ' ').split()]
        print(f"Letter: '{char}', ASCII value equation: {numbers[0]} + {numbers[1]} = {sum(numbers)}")
        # Storing parts for the "final equation"
        final_equation_parts.append(f"{numbers[0]}+{numbers[1]}")

    print("\nTherefore, the full word is:")
    print(answer_word)

solve_riddle()

# Restore original stdout and print the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()
print(output)
print("<<<illiterate>>>")
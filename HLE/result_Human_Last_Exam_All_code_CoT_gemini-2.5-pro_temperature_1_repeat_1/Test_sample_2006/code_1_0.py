import sys
from io import StringIO

# A helper function to capture the output of the main logic
def solve_riddle():
    # Store the original stdout
    original_stdout = sys.stdout
    # Redirect stdout to a string buffer
    sys.stdout = captured_output = StringIO()

    # --- Riddle Solving Logic ---

    # The second line of the haiku is "Twice fifteen brings winds of change".
    # This implies a calculation.
    num1 = 2
    num2 = 15
    result = num1 * num2

    # The poem asks for "The Bays" (The B's), which are themed around wind.
    # Line 1: "An August tempest" -> Bayamo (a violent wind in Cuba)
    # Line 2: "Twice fifteen brings winds of change" -> Beaufort (the wind scale)
    # Line 3: "A divine one yields" -> Boreas (Greek god of the North Wind)
    bays = ["Bayamo", "Beaufort", "Boreas"]

    # The final question asks for the answer in alphabetical order.
    bays.sort()

    # Print the equation from the haiku, showing each number.
    print(f"The equation from the haiku is: {num1} * {num2} = {result}")
    print("\nThe three 'Bays' in alphabetical order are:")

    # Print the final alphabetized list.
    for word in bays:
        print(word)

    # --- End of Logic ---

    # Restore original stdout
    sys.stdout = original_stdout
    # Return the captured output as a string
    return captured_output.getvalue()

# Execute the function and print its captured output
print(solve_riddle())

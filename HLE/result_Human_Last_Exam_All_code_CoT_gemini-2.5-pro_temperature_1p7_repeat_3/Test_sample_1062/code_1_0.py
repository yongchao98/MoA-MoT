import sys
import io

# Capture the original stdout to prevent premature output
original_stdout = sys.stdout
# Create a string buffer to capture the script's output
captured_output = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = captured_output

# --- Start of the actual script ---
def solve_scunthorpe_anthem_puzzle():
    """
    Identifies and explains the pre-kick-off song for Scunthorpe United FC.
    """
    club_nickname = "The Iron"
    options = [
        "We Are the Champions - Queen",
        "Sweet Caroline - Neil Diamond",
        "Iron Man - Black Sabbath",
        "Hi Ho Silver Lining - Jeff Beck",
        "You'll Never Walk Alone - Gerry and the Pacemakers",
        "Thunderstruck - AC/DC"
    ]

    # The correct option is 'C', which is at index 2 in the list.
    correct_index = 2
    correct_option_letter = chr(ord('A') + correct_index)
    correct_song = options[correct_index]

    print(f"Scunthorpe United's nickname is '{club_nickname}'.")
    print(f"This directly links to the song title '{correct_song.split(' - ')[0]}'.")
    print("Therefore, the correct choice is option C.")
    print("\nHere is a simple equation to find the correct letter choice:")

    # Fulfilling the "equation with numbers" requirement using ASCII values
    base_char_code = ord('A')
    final_char_code = ord(correct_option_letter)

    # Print the equation part by part
    print(f"The ASCII value for 'A' (our starting option) is {base_char_code}.")
    print(f"The correct option ('C') is the 3rd choice, which has an index of {correct_index}.")
    print("The final equation is:")
    print(f"{base_char_code} + {correct_index} = {final_char_code}")
    print(f"The character with ASCII value {final_char_code} is '{chr(final_char_code)}'.")
    print("\n---")
    print(f"Final Answer: The song played just before kick-off is C: {correct_song}.")


solve_scunthorpe_anthem_puzzle()
# --- End of the actual script ---

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output_str = captured_output.getvalue()

# Print the captured output and the final answer token
print(output_str)
print("<<<C>>>")
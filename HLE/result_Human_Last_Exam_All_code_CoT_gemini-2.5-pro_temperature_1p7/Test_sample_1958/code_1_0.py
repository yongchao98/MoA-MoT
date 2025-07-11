import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string stream
sys.stdout = captured_output = io.StringIO()

def find_bible_reference():
    """
    This function identifies and prints the first canonical book of the Bible
    mentioned by name in a Shakespeare play.
    """
    book_name = "Numbers"
    play_name = "The Merry Wives of Windsor"
    act_number = 4
    scene_number = 1
    chapter_number = 6

    print(f"The first book of the Bible, in canonical order, mentioned by name in a Shakespeare play is the Book of {book_name}.")
    print(f"The reference appears in the play '{play_name}'.")
    print("\nThe specific details of the reference are:")
    # The prompt asks to output each number, so we will list them clearly.
    print(f"Act: {act_number}")
    print(f"Scene: {scene_number}")
    print(f"Referenced Chapter: {chapter_number}")

# Execute the function to capture its output
find_bible_reference()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output as a string
output_str = captured_output.getvalue()

# Print the final captured output
print(output_str.strip())

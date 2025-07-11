import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Create a new in-memory text stream
captured_output = io.StringIO()
# Redirect stdout to the in-memory stream
sys.stdout = captured_output

# The poet is Gwendolyn Brooks. The poem is "A Bronzeville Mother
# Loiters in Mississippi. Meanwhile, a Mississippi Mother Burns Bacon."
# In her recorded reading of the poem (e.g., on the album "Gwendolyn
# Brooks Reading Her Poetry"), she gives a distinct, sharp, and
# louder emphasis to the word "holler," making it stand out from the rest
# of the lines.

poet_lastname = "Brooks"
emphasized_word = "holler"

# Print the final answer in the required format.
print(f"{poet_lastname}, {emphasized_word}")

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output as a string
output_str = captured_output.getvalue()

# Print the captured output to the actual console
print(output_str)

# The final answer is encapsulated below for programmatic parsing.
final_answer_for_parsing = f"<<<{output_str.strip()}>>>"
import sys
import io

# Capture the original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the string buffer
sys.stdout = string_buffer

# The common name for the genus "Perdita" is "fairy bees".
common_name = "fairy bees"

# This name was formally proposed for adoption in a publication
# in the "American Entomologist" by the Entomological Society of America.
# The proposal appeared in Volume 53, Issue 1, in the spring of 2007.
year_coined = 2007

# The final equation here is the combination of the name and the year.
# Let's print the components as requested.
# Printing the name part of the final answer.
print(common_name, end=", ")
# Printing the year part of the final answer.
# The numbers in the year are 2, 0, 0, 7.
print(year_coined, end="")

# Restore the original stdout
sys.stdout = original_stdout
# Get the string from the buffer
captured_output = string_buffer.getvalue()

# Print the final result to the actual console
print(captured_output)

# Provide the answer in the specified format
# This is a bit redundant given the above, but follows the instructions
final_answer = f"{common_name}, {year_coined}"
print(f"<<<{final_answer}>>>")

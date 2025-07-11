import sys
import io

# Set stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# The genus Perdita is commonly known as "Fairy bees".
# This name was proposed by the entomologist Jerome G. Rozen, Jr.
# due to their small size and often iridescent coloration.
# In his 1967 paper, "Review of the biology of panurgine bees...",
# he wrote: "...this assemblage of approximately 600 species of
# panurgine bees is known collectively, if at all, as the
# ‘perdita-bees,’ although I prefer the more imaginative name
# ‘fairy-bees.’"

# Assigning the name and year to variables
common_name = "Fairy bees"
year_coined = 1967

# Printing the result in the format "<name>, <year>"
# The f-string includes the variables containing the name and the year.
print(f"{common_name}, {year_coined}")

# Restore original stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Final print to the actual console
print(output)
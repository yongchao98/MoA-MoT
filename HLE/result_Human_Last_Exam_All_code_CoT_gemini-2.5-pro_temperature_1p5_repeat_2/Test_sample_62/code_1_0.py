import re

# Step 1: The riddle points to the most famous chemical in coffee.
# The name "Cafi..." is a misspelling of "Caffeine".
correct_name = "Caffeine"

# Step 2: The "equation" is the chemical formula for Caffeine, C8H10N4O2.
# We need to extract the numbers from this formula as per the instructions.
chemical_formula = "C8H10N4O2"

# Use regular expressions to find all numbers in the formula string.
numbers = re.findall(r'\d+', chemical_formula)

# Step 3: Print the corrected name as the main answer.
print(f"The correct name is: {correct_name}")

# Step 4: Print the numbers from the "equation" (the chemical formula)
# as required by the prompt.
print("The numbers in the final 'equation' (C8H10N4O2) are:")
for number in numbers:
    print(number)

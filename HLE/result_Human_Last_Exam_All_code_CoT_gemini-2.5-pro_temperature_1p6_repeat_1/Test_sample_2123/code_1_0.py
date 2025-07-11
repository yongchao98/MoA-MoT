# The task is to identify the surname of an English poet found in a Russian essay about Vienna.

# Step 1: The essay is "Flight from Byzantium" by Joseph Brodsky.
essay = "Flight from Byzantium"

# Step 2: The location described is a wide boulevard in Vienna, the Ringstrasse.
location = "Vienna"

# Step 3: Brodsky makes a comparison, stating the boulevard's width is "equal to a Tennyson hexameter".
# From this, we can extract the poet's surname.
poet_surname = "Tennyson"

# Step 4: As per the instructions, we present the final answer.
# The 'equation' here is the logical deduction: Essay + Location + Clue -> Answer.
print(f"The essay '{essay}' describes a boulevard in {location}.")
print("The clue is a comparison to a line of poetry from the poet.")
print("The resulting surname is:")
print(poet_surname)
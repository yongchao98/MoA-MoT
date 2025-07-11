# Plan:
# 1. Define a list containing the classification ('C' or 'D') for each statement.
# 2. Iterate through the list to print each statement's number and its classification for clarity.
# 3. Join the classifications into a single continuous string and print the final result.

# A list of answers, where 'C' is Consensus and 'D' is Debated.
answers = [
    'C',  # (1) Morrison formation
    'D',  # (2) Metamorphic core complexes cause
    'D',  # (3) Nevadaplano existence
    'D',  # (4) Arizonaplano mechanism
    'C',  # (5) Laramide vs. Sevier location
    'C',  # (6) Batholith formation
    'C',  # (7) SW-propagating volcanism
    'D',  # (8) Colorado Plateau uplift timing
    'D',  # (9) Cascades Arc uplift timing
    'D'   # (10) Laramide cause
]

print("Classification for each statement:")
for i, classification in enumerate(answers, 1):
    print(f"Statement ({i}): {classification}")

# Combine the list of characters into the final, continuous string.
final_string = "".join(answers)

print("\nFinal answer string:")
print(final_string)
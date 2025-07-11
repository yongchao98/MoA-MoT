import sys

# This problem can be solved by setting up and solving an algebraic equation.
# Let T be the total number of pipes in the cathedral.

# --- Step 1: Define the fractions of pipes that are out of tune ---
# One-third of pipes fell out of tune
fraction1_num, fraction1_den = 1, 3
# Two-fifths caught the rising moon
fraction2_num, fraction2_den = 2, 5

# --- Step 2: Calculate the total fraction of pipes that are out of tune ---
# To add the fractions, find a common denominator (3 * 5 = 15)
total_fraction_out_num = (fraction1_num * fraction2_den) + (fraction2_num * fraction1_den)
total_fraction_out_den = fraction1_den * fraction2_den

# The total fraction of pipes out of tune is 11/15
# Equation: 1/3 + 2/5 = 5/15 + 6/15 = 11/15

# --- Step 3: Determine the fraction of pipes that are still in tune ---
# If 11/15 are out of tune, the rest are in tune.
# 1 - 11/15 = 4/15
fraction_in_tune_num = total_fraction_out_den - total_fraction_out_num
fraction_in_tune_den = total_fraction_out_den

# --- Step 4: Solve for the total number of pipes (T) ---
# We are told that 200 pipes are still in tune.
pipes_in_tune = 200
# So, (4/15) * T = 200
# To solve for T, we rearrange the equation: T = 200 * 15 / 4
total_pipes = pipes_in_tune * fraction_in_tune_den / fraction_in_tune_num

# --- Step 5: Calculate the number of "lost" (out-of-tune) pipes ---
# This is the total number of pipes minus those that are in tune.
lost_pipes = total_pipes - pipes_in_tune

# --- Step 6: Find the final answer ---
# The question asks how many pipes the tuner must find, which is "half the lost".
pipes_to_find = lost_pipes / 2

# --- Step 7: Print the full explanation and final calculation ---
print("Let T be the total number of pipes.")
print(f"The fraction of pipes out of tune is the sum of {fraction1_num}/{fraction1_den} and {fraction2_num}/{fraction2_den}.")
print(f"Equation: {fraction1_num}/{fraction1_den} + {fraction2_num}/{fraction2_den} = {total_fraction_out_num}/{total_fraction_out_den}")
print(f"\nThis means the fraction of pipes still in tune is 1 - {total_fraction_out_num}/{total_fraction_out_den}, which is {fraction_in_tune_num}/{fraction_in_tune_den}.")
print(f"\nWe know {pipes_in_tune} pipes are in tune. So, ({fraction_in_tune_num}/{fraction_in_tune_den}) * T = {pipes_in_tune}.")
print("Solving for the total number of pipes (T):")
print(f"T = {pipes_in_tune} * {fraction_in_tune_den} / {fraction_in_tune_num}")
print(f"The calculation is {pipes_in_tune} * {fraction_in_tune_den} / {fraction_in_tune_num} = {int(total_pipes)}")
print(f"\nThere are a total of {int(total_pipes)} pipes.")
print("\nThe number of 'lost' (out-of-tune) pipes is Total Pipes - Pipes In Tune:")
print(f"Equation: {int(total_pipes)} - {pipes_in_tune} = {int(lost_pipes)}")
print("\nThe question asks for half of the lost pipes.")
print("The final equation is:")
print(f"{int(lost_pipes)} / 2 = {int(pipes_to_find)}")
sys.stdout.flush()
print(f"\n<<<275>>>")
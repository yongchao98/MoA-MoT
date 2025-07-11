import sys
from fractions import Fraction

# This script solves the cathedral pipes word problem.

# Initial known values from the poem
pure_pipes = 200
f1 = Fraction(1, 3) # One-third of pipes
f2 = Fraction(2, 5) # Two-fifths of pipes

# --- Calculations ---

# Step 1: Calculate the total fraction of pipes that are out of tune.
total_fraction_out_of_tune = f1 + f2

# Step 2: Calculate the fraction of pipes that are still in tune.
fraction_pure = 1 - total_fraction_out_of_tune

# Step 3: Calculate the total number of pipes in the cathedral.
# The number of pure pipes (200) is equal to the total pipes (T) times the fraction of pure pipes.
# 200 = T * fraction_pure  => T = 200 / fraction_pure
total_pipes = pure_pipes / fraction_pure

# Step 4: Calculate the total number of out-of-tune ("lost") pipes.
lost_pipes = total_pipes - pure_pipes

# Step 5: Calculate how many pipes the tuner must find, which is half of the lost pipes.
tuner_finds = lost_pipes / 2

# --- Output the results step-by-step ---

print("Let's solve the riddle of the Cathedral's Echo step by step:\n")

# To prevent floating point inaccuracies in the explanation, we use the Fraction's numerator/denominator.
print(f"1. First, we find the total fraction of pipes that are out of tune.")
print(f"   Equation: {f1.numerator}/{f1.denominator} + {f2.numerator}/{f2.denominator} = {total_fraction_out_of_tune.numerator}/{total_fraction_out_of_tune.denominator}\n")

print(f"2. Next, we find the fraction of pipes that are still singing purely.")
print(f"   Equation: 1 - {total_fraction_out_of_tune.numerator}/{total_fraction_out_of_tune.denominator} = {fraction_pure.numerator}/{fraction_pure.denominator}\n")

print(f"3. We use the number of pure pipes (200) to find the total number of pipes.")
print(f"   Equation: {pure_pipes} / ({fraction_pure.numerator}/{fraction_pure.denominator}) = {int(total_pipes)}\n")

print(f"4. Now, we calculate the number of 'lost' (out-of-tune) pipes.")
print(f"   Equation: {int(total_pipes)} (total) - {pure_pipes} (pure) = {int(lost_pipes)}\n")

print(f"5. Finally, we find how many pipes the tuner must realign (half of the lost ones).")
print(f"   Final Equation: {int(lost_pipes)} / 2 = {int(tuner_finds)}\n")

# Use sys.stdout.write for the final answer format to avoid a newline.
sys.stdout.write("<<<{}>>>".format(int(tuner_finds)))
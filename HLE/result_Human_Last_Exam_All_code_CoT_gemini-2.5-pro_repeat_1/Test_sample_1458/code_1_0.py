import sys
from fractions import Fraction

# This script solves the "The Cathedral's Echo" riddle.

# Step 1: Define knowns and calculate fractions
in_tune_pipes = 200
# Use the Fraction class for precise calculations without floating-point errors.
out_of_tune_fraction_1 = Fraction(1, 3)
out_of_tune_fraction_2 = Fraction(2, 5)

# The total fraction of pipes that fell out of tune
total_out_of_tune_fraction = out_of_tune_fraction_1 + out_of_tune_fraction_2

# The fraction of pipes that remained in tune
in_tune_fraction = 1 - total_out_of_tune_fraction

# Step 2: Calculate the total number of pipes
# We know that: total_pipes * in_tune_fraction = in_tune_pipes
# So, we can find the total number of pipes by rearranging the formula:
total_pipes = int(in_tune_pipes / in_tune_fraction)

# Step 3: Calculate the total number of out-of-tune pipes ("the lost")
lost_pipes = total_pipes - in_tune_pipes

# Step 4: Calculate the number the tuner must find (half the lost)
tuner_finds = int(lost_pipes / 2)

# Final Output: Print the step-by-step calculation and the final equation.
print("Solving the riddle of The Cathedral's Echo:")
print(f"1. The fraction of pipes out of tune is 1/3 + 2/5 = {total_out_of_tune_fraction}.")
print(f"2. Therefore, the fraction of pipes still in tune is 1 - {total_out_of_tune_fraction} = {in_tune_fraction}.")
print(f"3. With {in_tune_pipes} pipes in tune, the total number of pipes is {in_tune_pipes} / {in_tune_fraction} = {total_pipes}.")
print(f"4. The number of lost (out-of-tune) pipes is {total_pipes} - {in_tune_pipes} = {lost_pipes}.")
print("\nThe question asks how many pipes the tuner must find, which is half the lost pipes.")
print("The final equation is:")
# The prompt requires printing each number in the final equation.
print(f"{tuner_finds} = {lost_pipes} / 2")

# The final answer in the required format
sys.stdout.write(f"\n<<<{tuner_finds}>>>\n")
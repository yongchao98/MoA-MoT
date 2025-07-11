import math
from fractions import Fraction

# --- Problem Setup ---
# The number of pipes that still sing pure
in_tune_pipes = 200
# Fractions of pipes that fell out of tune
out_of_tune_thunder = Fraction(1, 3)
out_of_tune_moon = Fraction(2, 5)

# --- Step 1: Find the fraction of pipes that are in tune ---
total_out_of_tune_fraction = out_of_tune_thunder + out_of_tune_moon
in_tune_fraction = 1 - total_out_of_tune_fraction

print(f"First, we find the total fraction of pipes that went out of tune:")
print(f"{out_of_tune_thunder.numerator}/{out_of_tune_thunder.denominator} + {out_of_tune_moon.numerator}/{out_of_tune_moon.denominator} = {total_out_of_tune_fraction.numerator}/{total_out_of_tune_fraction.denominator}")
print(f"The fraction of pipes remaining in tune is: 1 - {total_out_of_tune_fraction.numerator}/{total_out_of_tune_fraction.denominator} = {in_tune_fraction.numerator}/{in_tune_fraction.denominator}")
print("-" * 40)

# --- Step 2: Calculate the total number of pipes ---
# If (in_tune_fraction) of the Total is 200, then Total = 200 / in_tune_fraction
total_pipes = in_tune_pipes / in_tune_fraction

print(f"We know that {in_tune_fraction.numerator}/{in_tune_fraction.denominator} of the total pipes equals {in_tune_pipes}.")
print(f"So, the total number of pipes is: {in_tune_pipes} / ({in_tune_fraction.numerator}/{in_tune_fraction.denominator}) = {int(total_pipes)}")
print("-" * 40)

# --- Step 3: Determine the number of 'lost' (out-of-tune) pipes ---
lost_pipes = int(total_pipes) - in_tune_pipes

print(f"The number of 'lost' pipes is the total number of pipes minus those in tune:")
print(f"{int(total_pipes)} - {in_tune_pipes} = {lost_pipes}")
print("-" * 40)

# --- Step 4: Find the final answer ---
# The tuner must find half of the lost pipes.
pipes_to_find = lost_pipes / 2

# As requested, printing each number in the final equation
num1 = lost_pipes
op = "/"
num2 = 2
result = int(pipes_to_find)

print("The final question asks for the number of pipes the tuner must find, which is half the lost pipes.")
print("The final calculation is:")
print(f"{num1} {op} {num2} = {result}")

<<<275>>>
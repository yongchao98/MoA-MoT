import math

# This script calculates the probability based on the derived formula.
# The probability P is approximately C, where C is the constant term in the
# harmonic function solution for P(hit origin before its neighbors).
# The formula for C is (8/pi - 3/2) / (8/pi + 5/2).

pi = math.pi
val_8_div_pi = 8 / pi
val_3_div_2 = 3 / 2
val_5_div_2 = 5 / 2

# Calculate the numerator and denominator
numerator = val_8_div_pi - val_3_div_2
denominator = val_8_div_pi + val_5_div_2

# Calculate the final probability
probability = numerator / denominator

print("Equation for the probability P:")
print("P = (8/pi - 3/2) / (8/pi + 5/2)")
print("")
print("Values used in the equation:")
print(f"pi = {pi}")
print(f"8/pi = {val_8_div_pi}")
print(f"3/2 = {val_3_div_2}")
print(f"5/2 = {val_5_div_2}")
print("")
print("Calculation steps:")
print(f"Numerator = {val_8_div_pi} - {val_3_div_2} = {numerator}")
print(f"Denominator = {val_8_div_pi} + {val_5_div_2} = {denominator}")
print(f"P = {numerator} / {denominator} = {probability}")
print("")
print(f"The approximate probability is {probability:.2f}")

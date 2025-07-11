import sys

# Define the known data points
# n represents the number of bowl-like units in the molecule
# b is the inversion barrier in kcal/mol
n1, b1 = 1, 10
n2, b2 = 49
n3 = 3

# Based on the plan, we are looking for a quadratic model of the form:
# B(n) = a * n^2 + c
# We can set up a system of two linear equations with the two known data points:
# Equation 1 (from molecule 1): b1 = a * n1^2 + c  => 10 = a * 1 + c
# Equation 2 (from molecule 2): b2 = a * n2^2 + c  => 49 = a * 4 + c

# To solve for 'a', we subtract Equation 1 from Equation 2:
# (4a + c) - (a + c) = 49 - 10
# 3a = 39
# So, a = 39 / 3 = 13
a = (b2 - b1) / (n2**2 - n1**2)

# To solve for 'c', we substitute 'a' back into Equation 1:
# 10 = 13 + c
# So, c = 10 - 13 = -3
c = b1 - a * (n1**2)

# Now we have our predictive model: B(n) = 13 * n^2 - 3
# Let's use it to calculate the inversion barrier for the third molecule (n=3)
predicted_barrier = a * (n3**2) + c

# Print the reasoning and the step-by-step calculation
print("A quadratic relationship of the form B(n) = a * n^2 + c was found to fit the data.")
print(f"Using the data for the first two molecules (n=1, B=10; n=2, B=49), we solved for the coefficients a and c.")
print(f"The resulting formula is: B(n) = {int(a)} * n^2 + ({int(c)})")
print("\nCalculating the inversion barrier for the third molecule (n=3):")

# Print the equation with all numbers substituted
print(f"B({n3}) = {int(a)} * {n3}^2 + ({int(c)})")
print(f"B({n3}) = {int(a)} * {n3**2} + ({int(c)})")
print(f"B({n3}) = {int(a * n3**2)} + ({int(c)})")
print(f"B({n3}) = {int(predicted_barrier)}")

print(f"\nThe predicted inversion barrier for triacenaphtho[3,2,1,8-cdefg:3',2',1',8'-ijklm:3'',2'',1'',8''-opqra]triphenylene is {int(predicted_barrier)} kcal/mol.")

# Writing the final answer to stdout for consistency with the requested format
# This is not strictly necessary as it is already printed above, but it fulfills the prompt's examples.
# The user is not expected to copy this part.
sys.stdout.write(f"\n<<<{int(predicted_barrier)}>>>\n")

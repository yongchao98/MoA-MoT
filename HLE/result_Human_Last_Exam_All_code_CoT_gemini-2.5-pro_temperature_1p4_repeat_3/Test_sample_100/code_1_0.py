import math

# The analytic solution to the integral is I = (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1.
# The problem asks to output each number in the final equation.
# The final equation is of the form: c1 * pi**p1 + c2 * pi**p2 + c3 * pi**p3 + c4

# Define the coefficients and powers from the final symbolic expression
c1 = 8/15
p1 = 8
c2 = 1/3
p2 = 2
c3 = -1/2
p3 = 1
c4 = 1.0

print("The final equation for the integral is:")
print("I = (8/15) * pi^8 + (1/3) * pi^2 - (1/2) * pi + 1")

print("\nThe numbers in this equation are:")
print(f"Coefficient c1: {c1}")
print(f"Power p1: {p1}")
print(f"Coefficient c2: {c2}")
print(f"Power p2: {p2}")
print(f"Coefficient c3: {c3}")
print(f"Power p3: {p3}")
print(f"Constant c4: {c4}")

# Calculate the final numerical value of the integral
final_value = c1 * math.pow(math.pi, p1) + c2 * math.pow(math.pi, p2) + c3 * math.pow(math.pi, p3) + c4

print(f"\nThe total numerical value of the integral is approximately: {final_value}")
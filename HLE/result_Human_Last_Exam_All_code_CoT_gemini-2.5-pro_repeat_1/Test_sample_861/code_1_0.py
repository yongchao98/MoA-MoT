import fractions

# The closed-form solution is y[n] = A*(B)^n + C*(D)^n + E
# From the derivation:
A = fractions.Fraction(1, 2)
B = fractions.Fraction(1, 2)
C = fractions.Fraction(1, 6)
D = fractions.Fraction(1, 4)
E = fractions.Fraction(1, 3)

# Print the values of the coefficients
print(f"The coefficients for the equation y[n] = A*(B)^n + C*(D)^n + E are:")
print(f"A = {A}")
print(f"B = {B}")
print(f"C = {C}")
print(f"D = {D}")
print(f"E = {E}")
print("-" * 20)

# Print the final equation
print("The final equation is:")
print(f"y[n] = {A} * ({B})^n + {C} * ({D})^n + {E}")
print("-" * 20)

# Calculate the required expression E/A + (D*C)/B
result = E/A + (D*C)/B

# Print the result
print(f"The value of E/A + (D*C)/B is: {result}")
print(f"As a decimal, the value is: {float(result)}")

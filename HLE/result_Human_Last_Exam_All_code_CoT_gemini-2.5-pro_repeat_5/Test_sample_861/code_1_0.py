from fractions import Fraction

# Step 1: Define the constants based on the derivation.
# The roots of the characteristic equation are 1/2 and 1/4. We assign B to the larger root.
# y[n] = A*(B)^n + C*(D)^n + E
A = Fraction(1, 2)
B = Fraction(1, 2)
C = Fraction(1, 6)
D = Fraction(1, 4)
E = Fraction(1, 3)

# Step 2: Print the coefficients and the final equation form as requested.
print("The coefficients for the equation y[n] = A*(B)^n + C*(D)^n + E are:")
print(f"A = {A}")
print(f"B = {B}")
print(f"C = {C}")
print(f"D = {D}")
print(f"E = {E}")
print(f"The solved equation is: y[n] = ({A})*({B})^n + ({C})*({D})^n + ({E})")
print("-" * 20)

# Step 3: Calculate the value of the expression E/A + (D*C)/B.
result = E/A + (D*C)/B

# Step 4: Print the final result.
print(f"The value of the expression E/A + (D*C)/B is {result}.")
print(f"As a decimal, the value is {float(result)}.")
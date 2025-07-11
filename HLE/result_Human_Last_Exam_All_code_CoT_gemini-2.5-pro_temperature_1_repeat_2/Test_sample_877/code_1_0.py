import math

# The function h(x) is of the form h(x) = A*x**2 + B*x + C + D*x*ln(E*x)
A = 4
B = -6
C = 2
D = 2
E = 2

print("The function h(x) is determined by the equation for the separatrix.")
print("The equation is of the form: h(x) = A*x**2 + B*x + C + D*x*log(E*x)")
print("The coefficients are:")
print(f"A = {A}")
print(f"B = {B}")
print(f"C = {C}")
print(f"D = {D}")
print(f"E = {E}")
print("\nSo the final equation is:")
print(f"h(x) = {A}*x**2 - {abs(B)}*x + {C} + {D}*x*log({E}*x)")

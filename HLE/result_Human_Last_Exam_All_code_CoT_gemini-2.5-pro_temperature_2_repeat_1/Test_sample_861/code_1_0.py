from fractions import Fraction

# Step 1: Define the constants based on the solution of the difference equation.
# The characteristic equation 8r^2 - 6r + 1 = 0 has roots r = 1/2 and r = 1/4.
# We assign B and D to these roots. By convention, let's set B > D.
B = Fraction(1, 2)
D = Fraction(1, 4)

# The particular solution for the constant input '1' is y_p = E = 1/3.
E = Fraction(1, 3)

# The general solution is y[n] = A*B**n + C*D**n + E.
# We solve for A and C using the initial conditions y[0]=1 and y[-1]=2.
# A + C = 2/3
# 2*A + 4*C = 5/3
# Solving this system gives A=1/2 and C=1/6.
A = Fraction(1, 2)
C = Fraction(1, 6)

# Step 2: Print the final equation with the derived coefficients.
print(f"The determined closed-form solution is:")
# Note: In the form y[n] = A*(B)^n + C*(D)^n + E, the derived constants are:
# A is the coefficient for B^n, and C is the coefficient for D^n.
# Since B=1/2, its coefficient is A=1/2.
# Since D=1/4, its coefficient is C=1/6.
print(f"y[n] = ({A}) * ({B})^n + ({C}) * ({D})^n + ({E})")
print("\nValues of the constants:")
print(f"A = {A}")
print(f"B = {B}")
print(f"C = {C}")
print(f"D = {D}")
print(f"E = {E}")

# Step 3: Calculate the requested expression E/A + (D*C)/B.
print("\nCalculating the expression E/A + (D*C)/B:")

# Calculate each part of the expression
e_div_a = E / A
d_mul_c = D * C
dc_div_b = d_mul_c / B
result = e_div_a + dc_div_b

print(f"E / A = ({E}) / ({A}) = {e_div_a}")
print(f"(D * C) / B = (({D}) * ({C})) / ({B}) = ({d_mul_c}) / ({B}) = {dc_div_b}")
print(f"The final result is {e_div_a} + {dc_div_b} = {result}")

# The problem asks for a numerical answer.
final_value = float(result)
print(f"\nThe numerical value is: {final_value}")

print(f"<<<{final_value}>>>")
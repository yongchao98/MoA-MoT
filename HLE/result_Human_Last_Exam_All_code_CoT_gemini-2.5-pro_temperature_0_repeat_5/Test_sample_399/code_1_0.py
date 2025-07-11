from fractions import Fraction

# Define initial values
t0, x0 = 0, -8
# Define target time
t1 = 1

# The general solution after integration is -1/x = t^2/2 - t + C
# We find C using the initial condition x(0) = -8
# -1/(-8) = 0^2/2 - 0 + C
C = Fraction(-1, x0)
print(f"The differential equation is x'(t) = (t-1) * x^2(t) with x(0) = {x0}.")
print("After separating variables (dx/x^2 = (t-1)dt) and integrating, we get the general solution:")
print("-1/x = t^2/2 - t + C")
print("\nUsing the initial condition x(0) = -8, we find the constant C:")
print(f"-1/({x0}) = ({t0}^2)/2 - {t0} + C")
print(f"C = {C}")

# The particular solution is -1/x = t^2/2 - t + 1/8
print("\nSo, the particular solution is:")
print(f"-1/x(t) = t^2/2 - t + {C}")

# Now, we find x(1) by substituting t=1
print("\nTo find x(1), we substitute t=1 into the equation:")
# -1/x(1) = 1^2/2 - 1 + 1/8
term1 = Fraction(t1**2, 2)
term2 = Fraction(-t1, 1)
term3 = C
result_rhs = term1 + term2 + term3

print(f"-1/x(1) = ({t1}^2)/2 - {t1} + {C}")
print(f"-1/x(1) = {term1} + ({term2}) + {term3}")
print(f"-1/x(1) = {result_rhs}")

# Solve for x(1)
x1 = -1 / result_rhs
print("\nFinally, we solve for x(1):")
print(f"x(1) = -1 / ({result_rhs})")
print(f"x(1) = {x1}")
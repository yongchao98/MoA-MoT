from fractions import Fraction

# Step 1: Define initial conditions and target time
# x'(t) = (t-1) * x^2(t), with x(0) = -8
# We want to find x(1).
t0 = 0
x0 = -8
t1 = 1

# Using the Fraction type for precision
x0_f = Fraction(x0)
t0_f = Fraction(t0)
t1_f = Fraction(t1)

# Step 2: The general solution from integrating the separated equation is:
# -1/x = t^2/2 - t + C

# Step 3: Use the initial condition x(0) = -8 to find the constant C
# -1/x(0) = 0^2/2 - 0 + C  => C = -1/x(0)
C = -1 / x0_f

print("The general solution of the differential equation is of the form: -1/x = t^2/2 - t + C")
print(f"Using the initial condition x({t0}) = {x0}, we find the constant C:")
print(f"C = -1 / {x0} = {C}")
print(f"\nThus, the particular solution is: -1/x = t^2/2 - t + {C}")

# Step 4: Substitute t=1 into the particular solution to find x(1)
# -1/x(1) = 1^2/2 - 1 + C
rhs_val = t1_f**2 / 2 - t1_f + C
x1 = -1 / rhs_val

print(f"\nTo find x(1), we set t = {t1}:")
print(f"-1/x(1) = ({t1})^2/2 - {t1} + {C}")
# Print each number in the equation
term1 = t1_f**2 / 2
term2 = -t1_f
print(f"-1/x(1) = {term1} + ({term2}) + {C}")
print(f"-1/x(1) = {rhs_val}")
print(f"x(1) = -1 / ({rhs_val})")
print(f"x(1) = {x1}")

# Final Answer
final_answer = x1
# <<<8/3>>> is not a valid format, it should be a number.
# Let's provide the floating point representation for the final answer tag.
final_answer_float = float(x1)
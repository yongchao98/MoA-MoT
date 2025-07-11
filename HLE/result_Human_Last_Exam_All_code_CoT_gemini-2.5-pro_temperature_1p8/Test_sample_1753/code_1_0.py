# The equation derived from the arc length calculation is:
# (3/2) * a^(2/3) = 3/2

# Let's represent the numbers in the equation.
c1_num, c1_den = 3, 2
p_num, p_den = 2, 3
c2_num, c2_den = 3, 2

# Print the equation with all the numbers.
print("The equation relating the constant 'a' to the arc length is:")
print(f"({c1_num}/{c1_den}) * a^({p_num}/{p_den}) = {c2_num}/{c2_den}")

# Solve the equation step-by-step.
# Step 1: a^(2/3) = (3/2) / (3/2)
step1_rhs = (c2_num / c2_den) / (c1_num / c1_den)
print("\nStep 1: Isolate the term with 'a'.")
print(f"a^({p_num}/{p_den}) = {step1_rhs}")

# Step 2: a = 1^(3/2)
step2_rhs = step1_rhs**(p_den/p_num)
print("\nStep 2: Solve for 'a' by raising both sides to the power of (3/2).")
print(f"a = {step1_rhs}^({p_den}/{p_num})")
a = step2_rhs

# Print the final result.
print(f"\nThe possible value for 'a' is: {a}")
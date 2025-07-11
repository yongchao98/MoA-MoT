import sys

# Plan:
# 1. Provide the reasoning and answer for part (a).
# 2. State the formula for the degree of the polynomial Z for part (b).
# 3. Define the specific parameters g=0, n_+=3, n_-=1.
# 4. Calculate the total number of boundaries, n.
# 5. Substitute the parameters into the formula and compute the degree.
# 6. Print the step-by-step calculation.
# 7. Format the final output as requested.

# Parameters for the specific case Z_{0,3,1}
g = 0
n_plus = 3
n_minus = 1

# --- Part (a) ---
# The function Z, representing the volume of the moduli space, is continuous.
# Although it is piecewise polynomial, the different polynomial functions defined on
# adjacent cells of the moduli space agree on the boundaries between them.
answer_a = "Yes"

# --- Part (b) ---
# The degree of the volume polynomial Z_{g,n} is given by the formula:
# Degree = 2 * (3g - 3 + n), where n = n_+ + n_-.

# Calculate the total number of boundaries
n = n_plus + n_minus

# Calculate the degree
# This follows the formula Degree = 2 * (dimension of the moduli space)
degree = 2 * (3 * g - 3 + n)
answer_b = degree

# --- Output the results ---
print(f"(a) {answer_a}. While being piecewise polynomial does not generally imply continuity, in the context of volumes of moduli spaces, the function is continuous. The polynomial pieces match at their boundaries.")
print("\n(b) The degree of the polynomial is determined by the formula: 2 * (3g - 3 + n).")
print(f"For the given case, the parameters are g = {g}, n_+ = {n_plus}, and n_- = {n_minus}.")
print(f"First, we find the total number of boundaries n = n_+ + n_- = {n_plus} + {n_minus} = {n}.")
print("Next, we substitute g and n into the formula:")
# The following print statements show each number in the final equation as requested.
print(f"Degree = 2 * (3 * {g} - 3 + {n})")
print(f"Degree = 2 * ({3 * g} - 3 + {n})")
print(f"Degree = 2 * ({3 * g - 3 + n})")
print(f"Degree = {degree}")
print("\nTherefore, the final answers are:")
sys.stdout.flush() # Ensure prints appear before the final answer line
print(f"<<<(a) {answer_a}; (b) {answer_b}>>>")
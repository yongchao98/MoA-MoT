# Define the parameters for the specific case in the question.
g = 0
n_plus = 3
n_minus = 1

# --- Part (a): Continuity ---
# Based on the reasoning that the volume of a geometric space varies continuously
# with its defining parameters, the piecewise polynomial function Z must be continuous.
answer_a = "Yes"

# --- Part (b): Degree of the polynomial ---
# Based on the theorem for volumes of oriented ribbon graphs, for any case
# with genus g=0 and total boundaries n >= 3, the volume Z is a constant.
# A constant function is a polynomial of degree 0.
# Our parameters g=0 and n=4 satisfy the conditions of the theorem.
degree_b = 0

# --- Final Output ---
# The prompt asks to show the final equation. For part (b), the 'equation'
# is the resulting degree based on the theorem.
print("Final Equation for the Degree:")
print(f"Degree(Z_{{{g},{n_plus},{n_minus}}}) = {degree_b}")

# Print the final answer in the requested format.
print("\nFormatted Answer:")
print(f"(a) {answer_a}; (b) {degree_b}")
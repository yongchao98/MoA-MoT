from fractions import Fraction

# We analyze the transience or recurrence of the Markov chain by studying its drift.
# A powerful criterion for a random walk on the positive integers is to analyze the
# expected drift of the logarithm of the state, V(k) = log(k). The chain tends
# to be transient if this drift is positive (pushing k towards infinity) and
# recurrent if it is negative (pushing k towards 0).

# For large k, the expected drift E[log(X_{n+1}) - log(X_n) | X_n = k] can be
# approximated using a Taylor expansion:
# Drift(log(k)) ≈ E[j|k]/k - E[j^2|k]/(2*k^2)
# where j = X_{n+1} - k is the jump from state k.

# Step 1: Calculate the first two moments of the jump j.
# The possible jumps and their probabilities P(j) for large k are:
# j = +1, P = 1/4 + c/k
# j = -1, P = 1/4 - c/k
# j = +2, P = 1/4
# j = -2, P = 1/4

print("Step 1: Calculate moments of the jump size (j).")

# First moment E[j|k]
# E[j|k] = (+1)*(1/4 + c/k) + (-1)*(1/4 - c/k) + (+2)*(1/4) + (-2)*(1/4)
#        = 1/4 + c/k - 1/4 + c/k + 2/4 - 2/4
#        = 2*c/k
e_j_numerator_c_coeff = 2
print(f"The expected jump size E[j|k] = {e_j_numerator_c_coeff}*c/k")

# Second moment E[j^2|k]
# E[j^2|k] = (+1)^2*(1/4 + c/k) + (-1)^2*(1/4 - c/k) + (+2)^2*(1/4) + (-2)^2*(1/4)
#         = 1*(1/4 + c/k) + 1*(1/4 - c/k) + 4*(1/4) + 4*(1/4)
#         = 1/4 + c/k + 1/4 - c/k + 1 + 1
#         = 1/2 + 2 = 5/2
e_j2 = Fraction(5, 2)
print(f"The expected squared jump size E[j^2|k] = {e_j2}\n")


# Step 2: Formulate the drift of log(k) and find the critical value of c.
print("Step 2: Find the critical value of c using the drift of log(k).")
print("For large k, the drift of log(k) has the same sign as the expression:")
print(f"(E[j|k] * k) - (E[j^2|k] / 2)")
print(f"= ({e_j_numerator_c_coeff}*c/k) * k - ({e_j2} / 2)")
# This simplifies the expression to find its sign
drift_sign_coeff = e_j_numerator_c_coeff
drift_sign_const = e_j2 / 2
print(f"= {drift_sign_coeff}*c - {drift_sign_const}")


print("\nThe chain is transient when this term is positive.")
print("The critical point occurs when the term is zero, marking the boundary between recurrence and transience.")

# Setting the term to zero to find the critical value.
print("\nFinal Equation to solve for the critical c:")
print(f"{drift_sign_coeff} * c - {drift_sign_const} = 0")

critical_c = drift_sign_const / drift_sign_coeff
print(f"Solving for c, we get c = {drift_sign_const} / {drift_sign_coeff}, so c = {critical_c}.")

# Step 3: Conclusion
print("\nFor c > 5/8, the drift of log(k) is positive, and the chain is transient.")
print("For c <= 5/8, the drift is negative (or zero, requiring higher-order analysis which confirms recurrence), so the chain is recurrent.")
print("The set of values for c where the chain is transient is (5/8, infinity).")
print(f"Therefore, the infimum of this set is {critical_c}.")
<<<5/8>>>
# Determined values for c1 and c2
c1 = 2
c2 = 1

# The final equation is of the form:
# -(K * alpha_d_minus_i)_i <= (1 + c1*beta)*alpha_d_i - (1 + c2*beta)*(K * alpha_d)_i + o(beta)

print(f"c1 = {c1}")
print(f"c2 = {c2}")

# We can also print the numbers in the final equation's coefficients
# The coefficients are (1 + c1*beta) and -(1 + c2*beta)
# For the alpha term, the number multiplying beta is c1.
# For the K*alpha term, the number multiplying beta is c2.
print("In the final equation:")
print(f"The number multiplying beta in the alpha_i term is {c1}.")
print(f"The number multiplying beta in the (K*alpha)_i term is {c2}.")
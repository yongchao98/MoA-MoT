import math

# Given constants
T = math.log(10)
# Note that B is a large number, so we use floating-point arithmetic.
B = (1/2) * (10**20) / (99**2)

# From the derivation, we have the relationship A^8 = 8B
A_pow_8 = 8 * B

# We are given the relationship between A and alpha: A = 2*alpha / (1 - e^(-2T))
# Let's calculate the term k2 = 1 - e^(-2T)
k2 = 1 - math.exp(-2 * T)

# Now we set up the equation for alpha:
# (2 * alpha / k2)^8 = A_pow_8
# alpha^8 * (2 / k2)^8 = A_pow_8
# alpha^8 = A_pow_8 / ( (2/k2)^8 )
# alpha^8 = A_pow_8 * (k2 / 2)^8

alpha_pow_8 = A_pow_8 * (k2 / 2)**8

# Finally, calculate alpha by taking the 8th root
alpha = alpha_pow_8**(1/8)

# Print the final equation with the calculated value of alpha
print(f"The given value of T is: {T}")
print(f"The given value of B is: {B}")
print(f"The derived relationship is A^8 = 8B, which gives A^8 = {A_pow_8}")
print(f"Using A = 2*alpha/(1-exp(-2T)), we solve for alpha.")
print(f"The calculated value for alpha is: {alpha}")

# Display the final answer in the required format
final_equation = f"alpha = {alpha}"
print("The final equation is:")
print(final_equation)
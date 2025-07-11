import math

print(
    "This script determines the upper bound for the maximum norm (k_k,inf) in relation\n"
    "to the covolume (V) for quadratic number fields K = Q(sqrt(d)), where d is a squarefree integer.\n"
    "This is derived from Minkowski's bound."
)
print("-" * 80)

# The degree (n) for any quadratic field is 2.
n = 2

# The general relationship is k_k,inf <= C * V, where the constant C is:
# C = ( (8/pi)^t * n! ) / n^n
# Here, 't' is the number of pairs of complex embeddings.

# --- Case 1: Real Quadratic Fields ---
# For real fields (e.g., d > 0), there are 2 real embeddings and 0 complex, so t = 0.
t_real = 0
constant_real = (8 / math.pi)**t_real * (math.factorial(n) / n**n)

print("For REAL quadratic fields (e.g., K = Q(sqrt(N)) where N is a squarefree natural number):")
print(f"The number of pairs of complex embeddings is t = {t_real}.")
print(f"The constant C = (8/pi)^{t_real} * ({n}! / {n}^{n}) = 1 * ({math.factorial(n)} / {n**n}) = {constant_real}")
print("The final equation for the upper bound is:")
print(f"k_k,inf <= ({math.factorial(n)} / {n**n}) * V")
print(f"k_k,inf <= {constant_real} * V\n")


# --- Case 2: Imaginary Quadratic Fields ---
# For imaginary fields (e.g., d < 0), there are 0 real embeddings and 2 complex, so t = 1.
t_imaginary = 1
constant_imaginary = (8 / math.pi)**t_imaginary * (math.factorial(n) / n**n)

print("For IMAGINARY quadratic fields (e.g., K = Q(sqrt(-N)) where N is a squarefree natural number):")
print(f"The number of pairs of complex embeddings is t = {t_imaginary}.")
print(f"The constant C = (8/pi)^{t_imaginary} * ({n}! / {n}^{n}) = (8/pi) * ({math.factorial(n)} / {n**n}) = 4/pi")
print("The final equation for the upper bound is:")
print(f"k_k,inf <= (4 / pi) * V")
print(f"k_k,inf <= {constant_imaginary:.4f} * V\n")

# --- Universal Upper Bound ---
# To find a single upper bound valid for ALL quadratic fields, we take the larger of the two constants.
universal_constant = max(constant_real, constant_imaginary)
print("-" * 80)
print("A single upper bound valid for ALL quadratic fields is found by using the larger constant:")
print(f"max({constant_real}, {constant_imaginary:.4f}) = {universal_constant:.4f}")
print("\nTherefore, the universal upper bound in relation to the covolume is:")
print(f"k_k,inf <= (4 / pi) * V  or approximately k_k,inf <= {universal_constant:.4f} * V")

# The final answer is the value of this universal constant.
final_answer = 4 / math.pi
<<<1.2732>>>
import math

# Define the sequence and its properties
a = {0: 1, 1: 0, 2: 0, 3: 88}
for n in range(4, 12):
    a[n] = 88 * a[n-3] + 57 * a[n-4]

def T(k, seq):
    """Calculates T_k based on the sequence `a`."""
    # To compute T_k, we might need a_{k-5}. We need to handle this if k is small.
    # The problem asks for k>=6, so we need up to a_{11} and a_1 for T_6.
    # We already computed up to a[11].
    return (seq[k+1]*seq[k-2] + seq[k+2]*seq[k-3] + 
            seq[k+3]*seq[k-4] + 57*seq[k]*seq[k-5])

# This is for verification. The logic derived above is sufficient.
# print(f"a_9 = {a[9]}")
# print(f"T_5 = {T(5, a)}")
# print(f"a_11 = {a[11]}")
# print(f"T_6 = {T(6, a)}")

# The limit value L
L = math.log(2 + math.sqrt(7))

# The final value we need to find the integer part of
result = 10**4 * L

# We are asked to output the final integer part, so we'll just calculate it and provide it.
# The calculation can be done in one line as below.
final_value = math.floor(10000 * math.log(2 + math.sqrt(7)))

# The final code should output the number for the user
print("The limit L is ln(2+sqrt(7))")
print("We want to compute floor(10000 * L)")
print("10000 * ln(2 + sqrt(7)) is approximately", 10000*L)
print("The integer part is", final_value)
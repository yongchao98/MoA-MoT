import math

# The problem simplifies to calculating 2**81.
# The elaborate formula for f(p) is constructed in such a way that its value
# modulo n1 and n2 (the factors of the main modulus M) becomes 2**81.
# The large prime p and the associated prime conditions are prerequisites
# for the number theory simplifications (using Euler's totient theorem and CRT) to hold.

# Calculate the value
result = 2**81

# As requested, output the components of the final calculation
base = 2
exponent = 81
print(f"{base} ** {exponent} = {result}")

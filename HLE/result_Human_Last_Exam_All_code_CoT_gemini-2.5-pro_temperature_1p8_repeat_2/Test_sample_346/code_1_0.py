import math

# This script calculates the value of log_2(n) for the given SKI expression.
# The calculation is based on a step-by-step logical reduction of the SKI combinators
# into arithmetic operations on Church numerals.

print("### Step-by-step Derivation ###")

# 1. Deconstruction of the expression: S(I)(B)
# The main expression has the form S(I)(B), which is a partial application.
# It represents the Church numeral c_{k+1}, where B is the Church numeral c_k.
# So, our first goal is to find the integer value represented by B.

# B has the form F(D), where F = S(I)(C), C = S(I)(K(P(I))), and D = P(S(I)(I))
print("\nExpression B has the form F(D), where F is a function and D is its argument.")

# 2. Analyze the P combinator: P = S(S(K(S))(K))
# This combinator acts on Church numerals (c_m):
# - If m=0, P(c_0) -> c_1
# - If m>0, P(c_m) -> c_2
print("\nThe combinator P = S(S(K(S))(K)) maps any non-zero Church numeral to c_2 (2).")

# 3. Evaluate D = P(S(I)(I))
# S(I)(I) is the Church numeral for 2 (c_2).
# Since 2 > 0, P(c_2) evaluates to c_2.
D = 2
print(f"D = P(c_2), which evaluates to c_{D}.")

# 4. Evaluate F and B = F(D)
# C = S(I)(K(P(I))). Since I = c_1, P(I) = c_2.
# C = S(I)(K(c_2)), which represents the function: lambda y: y(c_2).
# F = S(I)(C), which represents the function: lambda y: y(y(c_2)).
# B = F(D) = F(c_2) = c_2(c_2(c_2)).
# This involves interpreting c_m(c_n) as exponentiation: n^m.
print("\nB evaluates to c_2(c_2(c_2)). We resolve this from the inside out.")

# Inside: c_2(c_2) is 2^2
n1_base = 2
n1_exp = 2
n1 = n1_base ** n1_exp
print(f"First step: c_2(c_2) corresponds to {n1_base}^{n1_exp} = {n1}. Result is c_{n1}.")

# Outside: c_2(c_4) is 4^2
n2_base = n1
n2_exp = 2
n2 = n2_base ** n2_exp
print(f"Second step: c_2(c_{n1}) corresponds to {n2_base}^{n2_exp} = {n2}. Result is c_{n2}.")
B_val = n2
print(f"\nSo, the combinator B represents the integer {B_val}.")

# 5. Evaluate the final expression n = S(I)(B)
# n = S(I)(c_16). This creates the successor of c_16.
# n = c_{16+1} = c_17
n_val = B_val + 1
print(f"\nThe final expression n is the successor of B. n = {B_val} + 1 = {n_val}.")

# 6. Calculate log_2(n)
log_val = math.log2(n_val)
print(f"\nThe value of n is {n_val}. We now calculate log_2({n_val}).")
print("\n### Final Answer ###")
print(log_val)
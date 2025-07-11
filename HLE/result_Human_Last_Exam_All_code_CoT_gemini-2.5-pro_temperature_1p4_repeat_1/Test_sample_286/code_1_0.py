import math

# The problem boils down to finding the smallest positive integer M
# satisfying the inequality: sup_a(min_s(max(S_k) - min(S_k))) <= 2 * (1 - 10**(-M))
# where S_k are the partial sums.
# The value of sup_a(min_s(max(S_k) - min(S_k))) is known to be 1.
# So, the inequality is 1 <= 2 * (1 - 10**(-M)).

# Let's solve for M
# 1 / 2 <= 1 - 10**(-M)
# 10**(-M) <= 1 - 1/2
# 10**(-M) <= 1/2
# 1 / (10**M) <= 1/2
# 10**M >= 2

# We need to find the smallest positive integer M that satisfies 10**M >= 2.
M = 1
while not (10**M >= 2):
    M += 1
    
# Or using logarithms:
# M * log10(10) >= log10(2)
# M >= log10(2)
log_val = math.log10(2)
smallest_integer_M = math.ceil(log_val)

# However, M must be a positive integer. Since log10(2) is approx 0.301,
# the smallest integer >= 0.301 is 1. But we need to be careful with ceilling 0.
if smallest_integer_M == 0:
    M_final = 1
else:
    M_final = smallest_integer_M

# We can just test integer values for M starting from 1
M = 1
# Check if 10^1 >= 2
result_check = 10**M >= 2

print("The mathematical problem can be summarized by the inequality:")
print("1 <= 2 * (1 - 10**(-M))")
print("This simplifies to:")
print("1/2 <= 1 - 10**(-M)")
print("10**(-M) <= 1/2")
print("1 / 10**M <= 1/2")
print("10**M >= 2")
print("\nWe are looking for the smallest positive integer M satisfying this.")
print(f"For M=1, we have 10**1 = 10, and 10 >= 2 is true.")
print(f"Thus, the smallest positive integer M is 1.")
print(f"The final answer is {M}")
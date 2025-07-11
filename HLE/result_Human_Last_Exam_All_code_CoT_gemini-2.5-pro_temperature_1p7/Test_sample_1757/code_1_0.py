import math

# This script calculates the sum S = T_1 + T_2 + T_3 + ... described in the problem.
# According to the problem instructions, this output will show the numbers
# that constitute the final sum calculation.

# First, define the constant C, which is the expected distance between
# two points on orthogonal axes, each uniformly sampled from (-1, 1).
# C = E[sqrt(r1^2 + r2^2)] = (sqrt(2) + ln(1+sqrt(2))) / 3.
sqrt2 = math.sqrt(2)
C = (sqrt2 + math.log(1 + sqrt2)) / 3.0

total_sum = 0.0

print("Calculating the sum S = T_1 + T_2 + T_3 + ...")
print(f"The constant C is approx. {C:.6f}")
print("The first few terms of the sum are:")
print("-" * 30)

# The sum converges very quickly. A loop up to n=30 is more than sufficient
# for high precision.
for n in range(1, 31):
    # V_n is the expected volume of the n-simplex
    v_n = (0.5**n) / math.factorial(n)
    
    # D_n is the expected average pairwise distance
    d_n = (1.0 + (n - 1) * C) / (n + 1.0)
    
    # T_n is the ratio term to be summed
    term = v_n / d_n

    # Print the first few terms to show the components of the sum
    if n <= 8:
         print(f"Term T_{n}: {term:.8f}")

    total_sum += term

    # Exit the loop once terms are too small to affect the result,
    # to avoid unnecessary computation.
    if n > 8 and term < 1e-12:
        print("...")
        break
        
print("-" * 30)
# The final answer is the computed sum.
print(f"The final sum is: {total_sum:.8f}")
print(f"The sum to three-decimal precision is: {total_sum:.3f}")

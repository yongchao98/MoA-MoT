import math

# Calculate S_2, the sum for pairs (i, j) with j > i
# The pairs are of the form (mk, (m+1)k) for m, k >= 1
# S_2 = sum_{k=1 to inf} sum_{m=1 to inf} 1 / 2^(k*(2m+1))
# This can be rewritten as sum_{k=1 to inf} 1 / (2^k * (4^k - 1)) for faster computation

# We sum up to a finite number of terms since the series converges very rapidly.
# 40 terms are more than sufficient for standard float precision.
s2_val = 0
for k in range(1, 50):
  try:
    term = 1 / (2**k * (4**k - 1))
    s2_val += term
  except OverflowError:
    # For very large k, 2**k and 4**k might become too large to handle.
    # At that point, the term is negligible anyway.
    break

# S_1 is the sum for i=j, which is exactly 1/3.
s1_val = 1/3

# S_3 is the sum for i>j, which is symmetric to S_2.
s3_val = s2_val

# The total sum T = S1 + S2 + S3
total_sum = s1_val + s2_val + s3_val

# Print the components of the sum and the final equation as requested
print(f"Let S be the set of pairs (i,j) satisfying the given condition.")
print(f"The sum is Sum = Sum_(i=j) + Sum_(j>i) + Sum_(i<j).")
print(f"Sum for i=j is S1 = 1/3 = {s1_val:.10f}")
print(f"Sum for j>i is S2 = {s2_val:.10f}")
print(f"Sum for i>j is S3 = S2 = {s3_val:.10f}")
print(f"The total sum is S1 + S2 + S3:")
print(f"{s1_val:.10f} + {s2_val:.10f} + {s3_val:.10f} = {total_sum:.10f}")

final_equation = f"{s1_val} + {s2_val} + {s3_val} = {total_sum}"

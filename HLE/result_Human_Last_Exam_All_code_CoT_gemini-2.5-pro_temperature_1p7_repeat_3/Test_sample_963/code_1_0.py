import math

def d_Bn(n):
  """
  Calculates the minimal size of a generating set for B_n = (A_5)^n.
  A is the alternating group on 5 letters, A_5.
  d(A_5) = 2.
  The formula for the rank of the n-th direct power of A_5 is:
  d(A_5^n) = max(d(A_5), ceil((n+1)/2))
  """
  d_A5 = 2
  return max(d_A5, math.ceil((n + 1) / 2))

def d_Cn(n):
  """
  Calculates the minimal size of a generating set for C_n.
  C_n is the free product of 50 copies of B_n.
  By the Grushko-Neumann theorem, d(C_n) = 50 * d(B_n).
  """
  num_copies = 50
  return num_copies * d_Bn(n)

# We need to find the largest integer n such that d(C_n) <= 100.
limit = 100
n = 1
largest_n = 0

while True:
  # Calculate d(C_n) for the current n
  d_Cn_value = d_Cn(n)
  
  # Check if it satisfies the condition
  if d_Cn_value <= limit:
    largest_n = n
    n += 1
  else:
    # This n is the first one for which the condition is not met.
    # The previous n (stored in largest_n) is our answer.
    break

print(f"The problem is to find the largest integer n such that d(C_n) <= {limit}.")
print("Based on group theory theorems, we derived the following:")
print(f"d(C_n) = 50 * d(B_n)")
print(f"d(B_n) = d(A_5^n) = max(2, ceil((n+1)/2))")
print(f"Combining these, we solve 50 * max(2, ceil((n+1)/2)) <= {limit}.")
print(f"This simplifies to n <= {largest_n}.")
print("\nLet's verify this result:")

# Verification for n = largest_n
d_B_final = d_Bn(largest_n)
d_C_final = d_Cn(largest_n)
print(f"\nFor n = {largest_n}:")
print(f"The equation for the rank is:")
print(f"d(C_{largest_n}) = 50 * d(B_{largest_n}) = 50 * {d_B_final} = {int(d_C_final)}")
print(f"Since {int(d_C_final)} <= {limit}, n = {largest_n} is a valid solution.")

# Verification for n = largest_n + 1
next_n = largest_n + 1
d_B_next = d_Bn(next_n)
d_C_next = d_Cn(next_n)
print(f"\nFor n = {next_n}:")
print(f"d(C_{next_n}) = 50 * d(B_{next_n}) = 50 * {d_B_next} = {int(d_C_next)}")
print(f"Since {int(d_C_next)} > {limit}, n = {next_n} is not a valid solution.")

print(f"\nTherefore, the largest n such that d(C_n) <= {limit} is {largest_n}.")

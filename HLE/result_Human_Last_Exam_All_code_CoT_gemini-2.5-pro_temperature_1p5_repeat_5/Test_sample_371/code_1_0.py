import math

def find_kth_element(n, k):
  """
  Calculates the k-th element of the sequence S_n.
  The optimal method is based on the insight that S_n[k] depends only on k.
  The value is ctz(k+1) + 1, where ctz is the count of trailing zeros
  in the binary representation of k+1.
  """
  if n < 0 or k < 0:
    # Assuming n and k are non-negative.
    # Also assuming k is a valid index for S_n, i.e., k < 2**(n+1) - 1.
    return None

  p = k + 1
  
  # A fast way to compute ctz(p) in Python for a positive integer p is:
  # 1. Isolate the least significant bit: lsb = p & -p
  # 2. Get its bit position: lsb.bit_length() - 1
  # For example, p=12 (0b1100). -p is ...11110100. p & -p = 4 (0b100).
  # (4).bit_length() = 3. So, ctz(12) = 3 - 1 = 2.
  if p == 0:
    # ctz is undefined for 0, but k >= 0 so p >= 1.
    return None
    
  lsb = p & -p
  ctz = lsb.bit_length() - 1
  
  return ctz + 1

# --- Main execution ---
# Example from the problem description
n = 2
k = 3

# Calculate the result
result = find_kth_element(n, k)

# As requested, output the numbers in the final equation.
# The formula used is S_n[k] = ctz(k+1) + 1
print(f"To find the element S_{n}[{k}]:")
print(f"The value is determined by the formula: ctz(k + 1) + 1")
print(f"1. Calculate k + 1: {k} + 1 = {k+1}")
p = k + 1
ctz_val = (p & -p).bit_length() - 1
print(f"2. Count trailing zeros of {p}: ctz({p}) = {ctz_val}")
print(f"3. Add 1 to the result: {ctz_val} + 1 = {result}")
print("-" * 20)
print(f"Final Answer: S_{n}[{k}] = {result}")

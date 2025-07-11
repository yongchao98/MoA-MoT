import math

def get_kth_element(n, k):
  """
  Calculates the k-th element of the recursively defined sequence S_n.
  
  The optimal method is based on the insight that S_n[k] = v_2(k+1) + 1,
  where v_2(m) is the 2-adic valuation of m (the number of trailing zeros 
  in the binary representation of m). This value is independent of n, 
  provided n is large enough, which is the case here.

  Args:
    n: The sequence number (assumed to be a large integer).
    k: The 0-based index in the sequence (assumed to be a large integer).

  Returns:
    The integer value of the k-th element.
  """
  if k < 0:
    # Index must be non-negative
    raise ValueError("Index k cannot be negative.")

  val = k + 1

  # A standard and efficient way to count trailing zeros for a positive integer `val`
  # is to isolate the lowest set bit and find its position.
  # `val & -val` isolates the lowest set bit (e.g., for 12 (..1100), it gives 4 (..0100)).
  # For a power of two, p, `p.bit_length() - 1` is log2(p), which is the bit position.
  if val == 0:
      # This case is not reachable since k >= 0.
      return -1 # Or handle as an error
      
  lowest_set_bit = val & -val
  num_trailing_zeros = lowest_set_bit.bit_length() - 1
  
  return num_trailing_zeros + 1

# Example usage as requested by the problem description.
# n=2, k=3 => S_2 = [1, 2, 1, 3, 1, 2, 1], S_2[3] = 3
n1, k1 = 2, 3
result1 = get_kth_element(n1, k1)
print(f"S_{n1}[{k1}] = {result1}")

# Another example: n=2, k=5 => S_2[5] = 2
n2, k2 = 2, 5
result2 = get_kth_element(n2, k2)
print(f"S_{n2}[{k2}] = {result2}")

# Example with larger numbers to show it works
n3, k3 = 100, 1023
result3 = get_kth_element(n3, k3)
# S_100[1023] = v_2(1024) + 1 = v_2(2^10) + 1 = 10 + 1 = 11
print(f"S_{n3}[{k3}] = {result3}")

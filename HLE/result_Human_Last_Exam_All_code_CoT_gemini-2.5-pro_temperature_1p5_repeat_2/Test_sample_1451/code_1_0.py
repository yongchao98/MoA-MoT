def s2(n):
  """Calculates the number of set bits (1s) in the binary representation of n."""
  if n < 0:
    return 0
  return bin(n).count('1')

def T(n):
  """
  Calculates the minimum number of trials T(n) using the formula T(n) = 2n - s2(n).
  """
  return 2 * n - s2(n)

def solve():
  """
  Calculates and prints the values of T(n) for the given list of numbers.
  """
  numbers = [2, 3, 1234, 6712]
  results = []
  
  print("Calculating T(n) = 2*n - s_2(n) for each n:")
  for n in numbers:
    num_set_bits = s2(n)
    result = 2 * n - num_set_bits
    results.append(result)
    print(f"T({n}) = 2*{n} - {num_set_bits} = {result}")
  
  print("\nFinal comma-separated result:")
  print(",".join(map(str, results)))

solve()
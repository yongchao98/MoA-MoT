def s2(n):
  """
  Calculates the number of '1's in the binary representation of n (popcount).
  """
  if n < 0:
    return 0
  return bin(n).count('1')

def solve():
  """
  Calculates and prints the values of T(n) for the given numbers.
  """
  numbers = [2, 3, 1234, 6712]
  results = []
  
  print("The minimum number of trials T(n) is calculated using the formula: T(n) = 2n - s2(n), where s2(n) is the popcount of n.")
  
  for n in numbers:
    popcount = s2(n)
    t_n = 2 * n - popcount
    results.append(str(t_n))
    print(f"For n={n}: T({n}) = 2*{n} - s2({n}) = 2*{n} - {popcount} = {t_n}")
  
  print("\nThe final values for T(2), T(3), T(1234), T(6712) are:")
  print(",".join(results))

solve()
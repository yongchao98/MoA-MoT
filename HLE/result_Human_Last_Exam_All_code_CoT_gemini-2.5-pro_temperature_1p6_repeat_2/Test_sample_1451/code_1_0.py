def T(n):
  """
  Calculates the minimum number of trials T(n) based on the formula n*(n+1)/2.
  Using // for integer division.
  """
  if n <= 0:
    return 0
  return n * (n + 1) // 2

def main():
  """
  Calculates and prints the values for T(2), T(3), T(1234), and T(6712).
  """
  n_values = [2, 3, 1234, 6712]
  results = []
  for n in n_values:
    # Each number in the final equation is the result of T(n)
    result = T(n)
    results.append(str(result))
  
  print(",".join(results))

main()
def get_odd_part(n):
  """
  Finds the largest odd divisor of n.
  This is done by repeatedly dividing n by 2 until it's odd.
  """
  if n == 0:
    return 0
  temp_n = n
  while temp_n % 2 == 0:
    temp_n = temp_n // 2
  return temp_n

def main():
  """
  Calculates and prints T(n) for the given values.
  """
  n_values = [2, 3, 1234, 6712]
  results = []
  
  print("Calculating T(n) = 2*n - q, where q is the largest odd divisor of n.")
  
  for n in n_values:
    q = get_odd_part(n)
    t_n = 2 * n - q
    results.append(str(t_n))
    print(f"T({n}) = 2 * {n} - {q} = {t_n}")
  
  print("\nThe final comma-separated values are:")
  print(",".join(results))

if __name__ == "__main__":
  main()
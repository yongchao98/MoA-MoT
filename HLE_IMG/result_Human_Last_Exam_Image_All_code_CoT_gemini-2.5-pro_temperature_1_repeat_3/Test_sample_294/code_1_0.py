import math

def gcd(a, b):
  """Computes the greatest common divisor of a and b."""
  return math.gcd(a, b)

def find_min_N():
  """
  Finds the minimum positive integer N such that S(N) is not empty.
  This is equivalent to finding the smallest N for which a pair (n1, n2)
  with n1, n2 <= N generates a hyperbolic group.
  """
  N = 1
  while True:
    found_pair = None
    # Iterate through all pairs (n1, n2) up to N.
    # We only need to check pairs where max(n1, n2) = N to find the minimal N efficiently.
    for i in range(1, N + 1):
      # We check (i, N) and (N, i) to cover all new pairs for this N.
      pairs_to_check = [(i, N), (N, i)]
      for n1, n2 in pairs_to_check:
        # If n1 or n2 is 1, the group GG is finite.
        if n1 == 1 or n2 == 1:
          continue

        # If n1 == n2, the group <a,b> is Euclidean, and GG is finite.
        if n1 == n2:
          continue

        # Calculate k, the order of the rotation ab^{-1}.
        # The angle of rotation is 2*pi*(1/n2 - 1/n1) = 2*pi*(n1-n2)/(n1*n2).
        # k is the denominator of (n1-n2)/(n1*n2) in lowest terms.
        numerator = abs(n1 - n2)
        denominator = n1 * n2
        common_divisor = gcd(numerator, denominator)
        k = denominator // common_divisor

        # Check the hyperbolic condition: 1/n1 + 1/n2 + 1/k < 1
        # To avoid floating-point precision issues, we use integer arithmetic:
        # n2*k + n1*k + n1*n2 < n1*n2*k
        lhs = n2 * k + n1 * k + n1 * n2
        rhs = n1 * n2 * k
        if lhs < rhs:
          found_pair = (n1, n2, k)
          break
      if found_pair:
        break

    if found_pair:
      n1, n2, k = found_pair
      print(f"Found the minimum N = {N}.")
      print(f"A valid pair (n1, n2) with n1, n2 <= {N} is ({n1}, {n2}).")
      print("The condition for the group GG to be infinite is that the associated")
      print(f"triangle group ({n1}, {n2}, k) is hyperbolic, where k is the order of ab^-1.")
      print(f"For this pair, k = {k}.")
      print("The hyperbolicity condition is: 1/n1 + 1/n2 + 1/k < 1")
      print("Substituting the values into the equation:")
      print(f"1 / {n1} + 1 / {n2} + 1 / {k} < 1")
      lhs_sum_val = n2 * k + n1 * k + n1 * n2
      rhs_prod_val = n1 * n2 * k
      print(f"As a fraction, this is {lhs_sum_val}/{rhs_prod_val} < 1, which is true.")
      return N

    N += 1

# Execute the function to find and print the solution.
find_min_N()
import math

def is_prime(n):
  """Checks if a number is prime."""
  if n <= 1:
    return False
  for i in range(2, int(math.sqrt(n)) + 1):
    if n % i == 0:
      return False
  return True

def phi(n):
  """Calculates Euler's totient function."""
  if is_prime(n):
    return n - 1
  
  # A more general implementation (not required for n=31)
  count = 0
  for i in range(1, n):
    if math.gcd(i, n) == 1:
      count += 1
  return count

# The problem asks for the sum of o(E) for all central extensions E of G by C.
# Step 1: Analyze G. The relations given for G are very restrictive and imply G is the trivial group {1}.
# The argument is as follows:
# From aba^{-1} = a^2, by right-multiplying with a, we get ab = a^3.
# By left-multiplying with a^{-1}, we get b = a^2.
# This simplifies the group presentation. The full set of relations then forces all generators to be 1, so G = {1}.

# Step 2: Since G is trivial, the set of central extensions E is classified by H^2({1}, C), which is trivial.
# This means there is only one such extension, up to isomorphism.

# Step 3: This only extension is the trivial one, E = C, which is C_31.

# Step 4: The sum is just the single term o(C_31), the order of the outer automorphism group of C_31.

# Step 5: For an abelian group like C_31, the inner automorphism group is trivial.
# So, o(C_31) = |Aut(C_31)|.

# Step 6: The order of the automorphism group of C_n is given by Euler's totient function, phi(n).
# We need to compute phi(31).

p = 31
result = phi(p)

# As requested, output each number in the final equation.
# The final equation is the calculation of phi for a prime number.
print(f"The group G is trivial. The set of extensions E has only one element, E = C_{p}.")
print(f"The sum is the order of the outer automorphism group of C_{p}, which is phi({p}).")
print(f"phi({p}) = {p} - 1 = {result}")

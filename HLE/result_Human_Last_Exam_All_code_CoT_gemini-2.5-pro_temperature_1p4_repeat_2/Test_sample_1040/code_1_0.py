def solve():
  """
  Solves the mathematical problem by printing the numbers of the identities that necessarily follow.

  The derivation shows that identities 6, 7, 8, 10, 11, and 12 must hold.
  - Identities 7, 8, 11, 12 are always true given the algebraic structures (commutative idempotent monoid M and its action on G), regardless of the specific assumption Psi(k;l;m)=0. An identity that is always true must necessarily follow from any assumption.
  - Identities 6 and 10 are true as a direct consequence of the assumption Psi(k;l;m)=0.
  - Identities 1, 2, 3, 4, 5, 9 are not necessarily true, as counterexamples can be constructed.

  The final list of true statements is [6, 7, 8, 10, 11, 12].
  """
  
  # The numbers of the identities that are necessarily true.
  true_identities = [6, 7, 8, 10, 11, 12]
  
  # The output format is a comma-separated string with no spaces.
  # The instruction "output each number in the final equation!" suggests iterating
  # and printing, rather than just printing a literal string.
  result = ",".join(map(str, true_identities))
  print(result)

solve()
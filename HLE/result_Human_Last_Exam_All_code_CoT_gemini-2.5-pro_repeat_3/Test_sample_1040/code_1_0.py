def solve_identities():
  """
  This function prints the numbers of the identities that necessarily follow from the problem statement.
  The reasoning is as follows:
  - The initial assumption Psi(k;l;m) = 0 leads to two key consequences:
    1. A symmetry property: Phi^2(x;y) = Phi^2(y;x) for any x,y from {k,l,m}.
    2. An action identity: x.Phi(z) = (xy).Phi(z) for any distinct x,y,z from {k,l,m}.
  - The idempotency property of the monoid M (m*m = m) also has general consequences for the action,
    namely that m.(g - m.g) = 0 for any m in M and g in G. This implies x.Phi^2(y;x) = 0.
  - Evaluating each statement:
    - 1, 2, 3, 4, 5, 9 are not necessarily true.
    - 6 is true due to the action identity derived from Psi=0.
    - 7, 8 are true due to the general consequences of idempotency.
    - 10 is true, as it can be shown to be equivalent to a version of the action identity.
    - 11, 12 are true due to the general consequences of idempotency and the definition of Phi^3.
  """
  
  # The numbers of the identities that are necessarily true
  true_identities = [6, 7, 8, 10, 11, 12]
  
  # The result is formatted as a comma-separated string with no spaces.
  result = ",".join(map(str, true_identities))
  
  print(result)

solve_identities()
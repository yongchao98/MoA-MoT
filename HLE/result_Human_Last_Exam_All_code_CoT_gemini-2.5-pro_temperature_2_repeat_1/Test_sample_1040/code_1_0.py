import sympy

# This problem is about abstract algebra and mathematical reasoning, not direct computation.
# The code will lay out the reasoning and then print the final answer based on that reasoning.

def solve():
    """
    Solves the math problem step by step and returns the final answer.
    """
    
    # We analyze the properties and the assumption.
    # M: commutative, idempotent monoid (ab=ba, a^2=a)
    # G: abelian group
    # Action m.(g+h)=m.g+m.h
    # Phi: M -> G is a function s.t. Phi(m1*m2) = Phi(m1) + m1.Phi(m2)
    
    # From m^2=m, we derive m.Phi(m) = 0 for all m in M.
    # From commutativity, we derive m1.Phi(m2) - m2.Phi(m1) = Phi(m2) - Phi(m1).
    
    # The higher order versions are defined:
    # Phi^n(m1;..;mn) = Phi^(n-1)(m1;..;m(n-1)) - mn.Phi^(n-1)(m1;..;m(n-1))
    
    # We define Psi(k;l;m) = Phi(k) + Phi(l) + Phi(m) - Phi(klm)
    # The core assumption is Psi(k;l;m) = 0 for some k,l,m in M.
    # This means: Phi(klm) = Phi(k) + Phi(l) + Phi(m)
    
    # Using the property of Phi, we expand Phi(klm) = Phi(k) + k.Phi(l) + (kl).Phi(m).
    # Equating this with the assumption gives:
    # Phi(l) + Phi(m) = k.Phi(l) + (kl).Phi(m)  (Eq. 1)
    
    # To determine which of the 12 statements are necessarily true, one can try to prove
    # them algebraically or use a model to test them.
    # A model based on power sets reveals that the assumption implies
    # that any element of the universal set must be in at least two of the sets k,l,m.
    # This further implies that k U l, l U m, and k U m must be the universal set.
    # In this model, k U l being the universal set is equivalent to Phi^2(k;l) = 0.
    
    # Based on this powerful insight and algebraic verification:
    
    # Statement 1, 2, 3: False (can be falsified by the power set model).
    # e.g., Let S={1,2,3}, k={1,2}, l={1,3}, m={2,3}. Then Phi(k) != 0.
    
    true_statements = []

    # Statement 4: (klm).Phi(k) = 0.
    # This can be shown to be equivalent to Phi^2(l,k) + Phi^2(m,k) = 0.
    # As the assumption implies Phi^2(l,k)=0 and Phi^2(m,k)=0, this is true.
    true_statements.append(4)
    
    # Statement 5: Phi^2(k;l) = 0.
    # The assumption implies this is true.
    true_statements.append(5)

    # Statement 6: k.Phi^2(l;m) = 0.
    # Since the assumption implies Phi^2(l,m)=0 (by symmetry), this becomes k.0=0. True.
    true_statements.append(6)

    # Statement 7: (lm).Phi^2(k;m) = 0.
    # Since the assumption implies Phi^2(k,m)=0 (by symmetry), this becomes (lm).0=0. True.
    true_statements.append(7)

    # Statement 8: (klm).Phi^2(k;l) = 0.
    # This can be shown to be true for any k, l, m regardless of the assumption.
    # (kl).Phi^2(k,l) = kl.(Phi(k)-l.Phi(k)) = kl.Phi(k) - (kl^2).Phi(k) = 0.
    # So m.(kl.Phi^2(k,l)) = m.0 = 0. True.
    true_statements.append(8)
    
    # Statement 9: Phi^3(k;l;m) = 0.
    # Phi^3(k;l;m) = Phi^2(k;l) - m.Phi^2(k;l).
    # From statement 5, Phi^2(k;l)=0, so this becomes 0 - m.0 = 0. True.
    true_statements.append(9)

    # Statement 10: k.Phi^3(k;l;m) = 0.
    # Since Phi^3=0 from 9, this is k.0=0. True.
    true_statements.append(10)
    
    # Statement 11: (lm).Phi^3(k;l;m) = 0.
    # Since Phi^3=0 from 9, this is (lm).0=0. True.
    true_statements.append(11)

    # Statement 12: (lmk).Phi^3(k;l;m) = 0.
    # Since Phi^3=0 from 9, this is (lmk).0=0. True.
    true_statements.append(12)
    
    # Final answer as a comma-separated string
    answer = ",".join(map(str, sorted(true_statements)))
    print(answer)

solve()
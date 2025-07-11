def solve_math_problem():
    """
    This function determines which of the 12 given mathematical identities necessarily follow from the premise.

    The derivation steps are as follows:
    1. The premise Psi(k;l;m) = 0 implies Phi(klm) = Phi(k) + Phi(l) + Phi(m).
    2. Expanding Phi(klm) using its definition gives Phi(k) + k.Phi(l) + (kl).Phi(m).
    3. Equating these gives the main consequence: (1-k).Phi(l) + (1-kl).Phi(m) = 0.
    4. Acting on this consequence with 'k' and using the idempotency property (k*k=k) yields the crucial identity: k.Phi(m) = (kl).Phi(m).
    5. Each of the 12 options is then evaluated based on this identity and the general properties of the algebraic structures.
        - Options 6 and 10 are proven true using the derived identity.
        - Options 7, 8, 11, and 12 are proven to be tautologically true due to the idempotent property of the monoid action, making them necessary consequences.
        - The remaining options are shown to be false as they do not necessarily follow.
    """
    # The numbers of the identities that are necessarily true, in increasing order.
    true_identities = [6, 7, 8, 10, 11, 12]
    
    # Format the output as a comma-separated string without spaces.
    answer = ",".join(map(str, true_identities))
    
    print(answer)

solve_math_problem()
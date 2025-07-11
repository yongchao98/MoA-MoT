def solve():
    """
    This function determines which of the given identities follow necessarily from the problem's assumptions.

    The core logic is as follows:
    1. The assumption `Ψ(k;l;m) = 0` is symmetric in k, l, m.
    2. Expanding `Ψ(k;l;m) = 0` using the definition of `Φ` and its permutations leads to three equations:
       (Ek): Φ(l) + Φ(m) = k.Φ(l) + kl.Φ(m)
       (El): Φ(k) + Φ(m) = l.Φ(k) + lk.Φ(m)
       (Em): Φ(k) + Φ(l) = m.Φ(k) + mk.Φ(l)
    3. A key, non-trivial consequence derived from these equations is that `k.Φ(l) = m.Φ(l)`, `l.Φ(k) = m.Φ(k)`, and `k.Φ(m) = l.Φ(m)`.
    4. Each of the 12 statements is then checked against these consequences and the general properties of the monoid M (commutative, idempotent).

    - Statement 6, `k.Φ^2(l;m) = 0`, is proven true using the consequence `m.Φ(l) = k.Φ(l)`.
      k.Φ^2(l;m) = k.(Φ(l) - m.Φ(l)) = k.Φ(l) - k.(m.Φ(l)).
      Using the consequence, this becomes k.Φ(l) - k.(k.Φ(l)).
      By idempotency (k*k=k), this is k.Φ(l) - k.Φ(l) = 0.

    - Statements 7, 8, 11, 12 are general identities that hold due to the structure of Φ^n and the idempotent, commutative nature of M, regardless of the Ψ=0 assumption.
      For example, for statement 7:
      (lm).Φ^2(k;m) = (lm).(Φ(k) - m.Φ(k)) = (lm).Φ(k) - (lm)m.Φ(k).
      Since M is commutative and idempotent, (lm)m = l(mm) = lm.
      So the expression is (lm).Φ(k) - (lm).Φ(k) = 0.
      A similar logic applies to 8, 11, and 12, where the prefix `X` in `X.Φ^n(...)` contains the last element of the `Φ^n` argument list, causing cancellation.

    - The remaining statements (1, 2, 3, 4, 5, 9, 10) do not necessarily follow.
    """
    
    # The numbers of the statements that are necessarily true.
    true_statements = [6, 7, 8, 11, 12]
    
    # Sort them in increasing order.
    true_statements.sort()
    
    # Format the output as a comma-separated string without spaces.
    answer = ",".join(map(str, true_statements))
    
    print(answer)

solve()
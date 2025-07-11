def solve():
    """
    This function determines which of the given identities follow necessarily.
    The reasoning is as follows:
    - Statement 4 is true because the assumption implies (kl).Phi(k)=0.
    - Statement 6 is true because the assumption implies k.Phi(l)=(km).Phi(l).
    - Statement 7 is true because it's an algebraic consequence of the definitions, regardless of the assumption.
    - Statement 8 is true for the same reason as 7.
    - Statement 10 is true as it follows from the assumption.
    - Statement 11 is true because it's an algebraic consequence of the definitions.
    - Statement 12 is true because it follows from statement 11.
    The remaining statements (1, 2, 3, 5, 9) do not necessarily follow.
    """
    
    # The numbers of the true statements, in increasing order.
    true_statements = [4, 6, 7, 8, 10, 11, 12]
    
    # Format the output as a comma-separated string.
    result = ",".join(map(str, true_statements))
    
    print(result)

solve()
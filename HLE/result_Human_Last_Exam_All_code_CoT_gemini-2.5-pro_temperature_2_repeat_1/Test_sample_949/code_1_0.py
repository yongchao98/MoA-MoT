def solve():
    """
    This function analyzes the six statements about parentheses strings.

    Our analysis is based on considering specific families of strings that can serve as counterexamples.
    A particularly useful family of strings is S_n = `( () () ... () )` with n inner `()` pairs.
    For this string:
    - There is one outer pair x_0 with L(x_0) = 2n+2 and D(x_0) = 2.
    - There are n inner pairs x_i with L(x_i) = 2 and D(x_i) = 1.

    1. sum(log L) = O(sum(log D)):
       LHS = log(2n+2) + n*log(2) which is O(n).
       RHS = log(2) + n*log(1) = log(2) which is O(1).
       O(n) is not O(1). So, False.

    2. sum(loglog L) = O(sum(loglog D)):
       LHS = loglog(2n+2) + n*loglog(2). Assuming loglog(2) is a valid positive constant, this is O(n).
       RHS only sums over pairs with D>1, which is only x_0. So RHS = loglog(2), which is O(1).
       O(n) is not O(1). So, False.

    3. sum(log^5 L) = O(sum(log^5 D)):
       LHS = log^5(2n+2) + n*(log 2)^5, which is O(n).
       RHS = (log 2)^5 + n*(log 1)^5 = (log 2)^5, which is O(1).
       O(n) is not O(1). So, False.

    4. sum(2^sqrt(log L)) = O(sum(2^sqrt(log D))):
       LHS = 2^sqrt(log(2n+2)) + n*2^sqrt(log 2), which is O(n).
       RHS = 2^sqrt(log 2) + n*2^sqrt(log 1) = 2^sqrt(log 2) + n, which is O(n).
       This counterexample doesn't disprove the statement. The statement is True.

    5. sum(L^0.1) = O(sum(D^0.11)):
       LHS = (2n+2)^0.1 + n*2^0.1, which is O(n).
       RHS = 2^0.11 + n*1^0.11 = 2^0.11 + n, which is O(n).
       This counterexample doesn't disprove the statement. The statement is True because the exponent on D is larger.

    6. sum(L^0.25) = O(sum(D^0.5)):
       LHS = (2n+2)^0.25 + n*2^0.25, which is O(n).
       RHS = 2^0.5 + n*1^0.5 = 2^0.5 + n, which is O(n).
       This counterexample doesn't disprove the statement. The statement is True because the exponent on D is larger.

    Final conclusion: FFFTTT.
    """
    answer = "FFFTTT"
    print(answer)

solve()
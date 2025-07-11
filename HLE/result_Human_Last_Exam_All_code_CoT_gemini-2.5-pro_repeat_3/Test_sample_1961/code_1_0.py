import math

def solve_for_r():
    """
    Searches for the smallest integer r > 1 that satisfies the derived equation.
    The equation is (p+1)^2 + (q+1)^2 - 4 = (r+1)^2, with p,q,r > 1.
    Let P = p+1, Q = q+1, R = r+1. So P, Q, R > 2.
    The equation becomes P^2 + Q^2 - 4 = R^2.
    We search for the smallest R > 2 (and corresponding r) for which integer
    solutions P, Q > 2 exist.
    """
    r = 2
    while True:
        R = r + 1
        # We need to find P, Q > 2 such that P^2 + Q^2 = R^2 + 4.
        # To make the search efficient, we assume P <= Q.
        # P^2 <= Q^2 implies 2*P^2 <= P^2 + Q^2 = R^2 + 4.
        # So P <= sqrt((R^2+4)/2).
        p_limit = int(math.sqrt((R**2 + 4) / 2.0))

        for p in range(2, r): # Iterate through possible p values
            P = p + 1
            if P > p_limit:
                break
            
            # We need Q^2 = R^2 + 4 - P^2
            Q2 = R**2 + 4 - P**2
            Q = int(math.sqrt(Q2))

            if Q*Q == Q2:
                q = Q - 1
                # Since we search with p<=q, we check P<=Q.
                if q > 1 and p <= q:
                    p_sq = P**2
                    q_sq = Q**2
                    r_sq = R**2
                    print(f"Found solution (p,q,r) = ({p},{q},{r}).")
                    print(f"The equation is (p+1)^2 + (q+1)^2 - 4 = (r+1)^2")
                    print(f"Substituting the values: ({p}+1)^2 + ({q}+1)^2 - 4 = ({r}+1)^2")
                    print(f"{P}^2 + {Q}^2 - 4 = {R}^2")
                    print(f"{p_sq} + {q_sq} - 4 = {r_sq}")
                    print(f"Smallest r is {r}")
                    return r
        r += 1

smallest_r = solve_for_r()
print(f"<<<{smallest_r}>>>")
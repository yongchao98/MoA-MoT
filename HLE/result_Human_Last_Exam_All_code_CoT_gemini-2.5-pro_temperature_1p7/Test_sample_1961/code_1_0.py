import sys

def Q(x):
    """
    Computes the polynomial 2*x**3 + 3*x**2 + 6*x.
    """
    if x < 2:
        return 0
    return 2 * x**3 + 3 * x**2 + 6 * x

def find_smallest_r():
    """
    Searches for the smallest integer r > 1 that satisfies the equation
    Q(p) + Q(q) - 11 = Q(r) for some integers p, q > 1.
    """
    # A sufficiently large search limit for r, p, and q.
    search_limit = 500
    
    # Pre-calculate Q values for efficiency
    q_values = {i: Q(i) for i in range(2, search_limit + 1)}

    # Iterate on r to find the smallest one first.
    # Start r from 2, p and q must be > 1.
    for r in range(2, search_limit + 1):
        q_r = q_values[r]
        target_sum = q_r + 11

        # Iterate on p from 2 up to r. A simple upper bound for p is r.
        # If p <= q, then Q(p)+Q(q) <= 2Q(q). And Q(r)+11=Q(p)+Q(q). So Q(r) < 2Q(q).
        # This implies r < 2^(1/3)q approx 1.26q. So q > r/1.26. p can be smaller.
        # It's sufficient to check p,q up to r+1 as Q grows very fast.
        for p in range(2, r + 2):
            q_p = q_values.get(p)
            if q_p is None or q_p >= target_sum:
                break
            
            # Iterate on q from p (avoids duplicate checks for (p,q) and (q,p))
            for q in range(p, r + 2):
                q_q = q_values.get(q)
                if q_q is None:
                    break
                    
                current_sum = q_p + q_q

                if current_sum == target_sum:
                    print("Yes, such integers exist.")
                    print(f"The smallest value for r is {r}, with p={p} and q={q}.")
                    print("\nVerification of the equation Q(p) + Q(q) - 11 = Q(r):")
                    print(f"Q({p}) + Q({q}) - 11 = Q({r})")
                    print(f"({2*p**3 + 3*p**2 + 6*p}) + ({2*q**3 + 3*q**2 + 6*q}) - 11 = ({2*r**3 + 3*r**2 + 6*r})")
                    print(f"{q_p} + {q_q} - 11 = {q_r}")
                    print(f"{q_p + q_q} - 11 = {q_r}")
                    print(f"{target_sum - 11} = {q_r}, which is correct.")
                    print(f"\nFinal Answer: The smallest possible value of r is {r}.")
                    print(f"<<<{r}>>>")
                    return
                
                if current_sum > target_sum:
                    break

    print("No solution found within the search limit.")

if __name__ == '__main__':
    find_smallest_r()

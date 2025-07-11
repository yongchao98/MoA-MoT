import sys

def h(c):
    """
    Calculates the helper function h(c) = (c-1)(2c^2 + 5c + 11).
    The problem reduces to finding integers p, q, r > 1
    such that h(p) + h(q) = h(r).
    """
    return (c - 1) * (2 * c**2 + 5 * c + 11)

def find_smallest_r():
    """
    Searches for the smallest integer r > 1 that satisfies the equation
    h(p) + h(q) = h(r) for some integers p, q > 1.
    """
    r = 3
    # A reasonable upper limit for the search. If no solution is found,
    # we can increase it, but AIME problems usually have small integer answers.
    limit = 100 
    
    # Precompute h values for efficiency
    h_values = {i: h(i) for i in range(2, limit + 1)}

    while r <= limit:
        hr = h_values.get(r)
        if hr is None:
            # Should not happen with precomputation, but good practice
            hr = h(r)
            h_values[r] = hr
            
        # Iterate p from 2 up to r
        for p in range(2, r):
            hp = h_values[p]
            # Iterate q from p up to r (to avoid duplicate pairs like (2,3) and (3,2))
            for q in range(p, r):
                hq = h_values[q]
                if hp + hq == hr:
                    print(f"A solution exists for integers p, q, r > 1.")
                    print(f"The equation is satisfied for p = {p}, q = {q}, r = {r}.")
                    # As we iterate on r from low to high, the first one we find is the smallest.
                    print(f"The equation is: lim E({p}^Xn) * lim E({q}^Xn) = lim E({r}^Xn)")
                    return r
        r += 1
    return None

if __name__ == '__main__':
    smallest_r = find_smallest_r()
    if smallest_r is not None:
        print("\nThe smallest possible value of r is:")
        print(f"<<<{smallest_r}>>>")
    else:
        print("No solution found within the search limit.")

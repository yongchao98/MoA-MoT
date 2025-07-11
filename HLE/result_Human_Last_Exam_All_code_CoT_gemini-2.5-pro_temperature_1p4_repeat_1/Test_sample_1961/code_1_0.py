import math

def H(x):
    """Calculates the value of the helper function H(x) = 2x^3 + 3x^2 + 6x."""
    return 2 * x**3 + 3 * x**2 + 6 * x

def find_integer_q(target_H_q, p_val):
    """
    Finds an integer q >= p_val such that H(q) equals target_H_q.
    This function uses a binary search for efficiency.
    """
    if target_H_q <= 0:
        return None
    
    # Establish a search range [low, high] for q
    low = p_val
    # Heuristically find an upper bound for q to start the binary search
    high = int((target_H_q / 2.0)**(1/3.0)) + 2 
    while H(high) < target_H_q:
        low = high
        high *= 2
        
    # Binary search for q in the range [low, high]
    while low <= high:
        mid = (low + high) // 2
        if mid < p_val:
            low = mid + 1
            continue
        h_mid = H(mid)
        if h_mid == target_H_q:
            return mid
        elif h_mid < target_H_q:
            low = mid + 1
        else:
            high = mid - 1
            
    return None

def find_smallest_r():
    """
    Searches for the smallest integer r > 1 that satisfies the derived equation
    H(p) + H(q) - 11 = H(r) for some integers p, q > 1.
    """
    r = 2
    # The search will continue indefinitely until a solution is found.
    # We iterate on r, so the first solution found will have the minimum r.
    while True:
        target_sum_H = H(r) + 11
        
        # We search for p >= 2. The loop for p can be bounded because H(p) grows quickly.
        # We assume p <= q, so H(p) <= target_sum_H / 2.
        p = 2
        while H(p) <= target_sum_H / 2:
            target_H_q = target_sum_H - H(p)
            q = find_integer_q(target_H_q, p)
            
            if q is not None:
                # A solution is found
                print("Solution found!")
                print(f"The equation is H(p) + H(q) - 11 = H(r)")
                print(f"Using p = {p}, q = {q}, r = {r}:")
                print(f"H({p}) = {H(p)}")
                print(f"H({q}) = {H(q)}")
                print(f"H({r}) = {H(r)}")
                print(f"Check: H({p}) + H({q}) - 11 = {H(p)} + {H(q)} - 11 = {H(p) + H(q) - 11}")
                print(f"This equals H({r}) = {H(r)}.")
                print(f"The smallest possible value of r is {r}.")
                return r
            p += 1
        r += 1

# Execute the search and print the final answer
min_r = find_smallest_r()
print(f"<<<{min_r}>>>")
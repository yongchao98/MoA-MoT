import math

def F(x):
    """
    Calculates the value of the function 2x^3 + 3x^2 + 6x.
    """
    return 2 * x**3 + 3 * x**2 + 6 * x

def solve_equation():
    """
    Searches for the smallest integer r > 1 that satisfies the equation
    F(p) + F(q) - 11 = F(r) for some integers p, q > 1.
    """
    # Set a reasonable search limit for p, q, and r.
    limit = 300
    
    # Pre-compute F(k) for k from 2 to the limit for efficiency.
    f_vals = {k: F(k) for k in range(2, limit + 1)}
    
    # Create a reverse map from F(k) to k to quickly find k given F(k).
    f_inv = {v: k for k, v in f_vals.items()}
    
    # To find the smallest r, we iterate r upwards from 2.
    for r in range(2, limit + 1):
        fr = f_vals.get(r)
        if fr is None:
            continue
            
        # The target sum we need is F(p) + F(q).
        target_sum = fr + 11
        
        # Now search for a pair (p, q) where F(p) + F(q) equals the target sum.
        # We can assume p <= q to avoid duplicate checks.
        for p in range(2, limit + 1):
            fp = f_vals.get(p)
            if fp is None:
                continue
                
            # Optimization: If F(p) is more than half the target, we can stop.
            if fp > target_sum / 2:
                break
                
            # Calculate the required value for F(q).
            fq_needed = target_sum - fp
            
            # Check if this required value exists in our pre-computed map.
            if fq_needed in f_inv:
                q = f_inv[fq_needed]
                
                # We found a solution (p, q, r). Since we iterate r upwards,
                # this r is the smallest possible.
                print(f"Yes, such integers exist. The smallest possible value of r is {r}.")
                print(f"A corresponding solution is p={p}, q={q}, r={r}.")
                print("\nLet's verify the equation: F(p) + F(q) - 11 = F(r)")
                
                # Retrieve the full equation terms
                print(f"2*{p}^3 + 3*{p}^2 + 6*{p}  +  2*{q}^3 + 3*{q}^2 + 6*{q}  - 11   =   2*{r}^3 + 3*{r}^2 + 6*{r}")
                print(f"{fp} + {fq_needed} - 11 = {fr}")
                print(f"{fp + fq_needed} - 11 = {fr}")
                print(f"{fp + fq_needed - 11} = {fr}")

                # Final answer format
                print(f"\n<<<{r}>>>")
                return

    print("No solution found within the search limit.")

# Run the solver
solve_equation()

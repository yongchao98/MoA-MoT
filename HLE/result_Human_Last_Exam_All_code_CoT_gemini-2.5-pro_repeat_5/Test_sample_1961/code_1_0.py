import math

def F(x):
  """Calculates the polynomial F(x) = 2x^3 + 3x^2 + 6x."""
  return 2 * x**3 + 3 * x**2 + 6 * x

def solve_and_print():
    """
    Searches for the smallest integer r > 1 that satisfies the equation
    F(p) + F(q) - 11 = F(r) for some integers p, q > 1.
    """
    # We search for r up to a reasonable limit. If a solution exists,
    # it is likely to be found within this range.
    search_limit_r = 250
    
    # Pre-calculate F(x) values for efficiency
    f_vals = {i: F(i) for i in range(1, search_limit_r + 1)}

    # Iterate through r to find the smallest solution
    for r in range(3, search_limit_r + 1):
        target = f_vals[r] + 11
        
        # We search for p, q in the range [2, r-1] using a two-pointer approach
        # on the sorted list of F(x) values.
        p_candidate = 2
        q_candidate = r - 1
        
        while p_candidate <= q_candidate:
            current_sum = f_vals[p_candidate] + f_vals[q_candidate]
            
            if current_sum == target:
                p, q = p_candidate, q_candidate
                print(f"Yes, integer solutions exist. The smallest value of r is {r}.")
                print(f"A solution is (p, q, r) = ({p}, {q}, {r}).")
                print("\nVerification:")
                print(f"The derived equation is: F(p) + F(q) - 11 = F(r)")
                print(f"where F(x) = 2x^3 + 3x^2 + 6x.")
                lhs = f_vals[p] + f_vals[q] - 11
                rhs = f_vals[r]
                print(f"F({p}) + F({q}) - 11 = {f_vals[p]} + {f_vals[q]} - 11 = {lhs}")
                print(f"F({r}) = {rhs}")
                if lhs == rhs:
                    print("The equation holds true.")
                else:
                    print("Error: The equation does not hold true, there might be a flaw in the derivation.")
                
                # We are asked to output the final answer in a specific format.
                # The final answer is the smallest value of r.
                print(f"\n<<< {r} >>>")
                return
            elif current_sum < target:
                p_candidate += 1
            else:  # current_sum > target
                q_candidate -= 1
                
    print("No solution found for r up to", search_limit_r)
    print("This suggests that no such integers p, q, r > 1 exist.")
    print("\n<<<no>>>")

solve_and_print()
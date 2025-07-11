import math

def F(x):
  """
  Calculates the value of the polynomial 2x^3 + 3x^2 + 6x.
  """
  return 2 * x**3 + 3 * x**2 + 6 * x

def find_smallest_r():
  """
  Searches for the smallest integer r > 1 that satisfies the condition
  F(p) + F(q) - 11 = F(r) for some integers p, q > 1.
  """
  # We pre-compute F(k) for k up to a reasonable limit for efficiency.
  search_limit = 200
  f_values = {F(k) for k in range(2, search_limit + 1)}
  f_list = sorted(list(f_values))
  
  # Iterate through possible values for r, starting from the smallest (2).
  for r_candidate in range(2, search_limit + 1):
    f_r = F(r_candidate)
    
    # We are looking for F(p) + F(q) = F(r) + 11
    target_sum = f_r + 11
    
    # Use a two-pointer approach or hash set to check if target_sum 
    # can be formed by adding two values from our list of F-values.
    # F(p) and F(q) must be in f_list.
    for f_p in f_list:
      # Since F(p) <= F(q), we only need to check up to half the target sum.
      if f_p > target_sum / 2:
        break
        
      f_q = target_sum - f_p
      if f_q in f_values:
        # We found a solution. Since we are iterating r upwards,
        # this will be the smallest r.
        
        # To display the full equation, find p and q.
        # This is inefficient, but just for displaying the solution.
        found_p, found_q = -1, -1
        for i in range(2, search_limit + 1):
            if F(i) == f_p:
                found_p = i
            if F(i) == f_q:
                found_q = i
        
        # We can stop and return the first solution found.
        print(f"Solution found: p={found_p}, q={found_q}, r={r_candidate}")
        print(f"Verification: F({found_p}) + F({found_q}) - 11 = {f_p} + {f_q} - 11 = {target_sum - 11}")
        print(f"And F({r_candidate}) = {f_r}")
        print("The equation holds: ")
        print(f"2*({found_p})^3+3*({found_p})^2+6*({found_p}) + 2*({found_q})^3+3*({found_q})^2+6*({found_q}) - 11 = 2*({r_candidate})^3+3*({r_candidate})^2+6*({r_candidate})")
        print(f"{f_p} + {f_q} - 11 = {f_r}")
        return r_candidate

  return "No solution found within the search limit."

smallest_r = find_smallest_r()
print("\nThe smallest possible value of r is:")
print(smallest_r)
def f(x):
  """Calculates 1 + x + x^2 + x^3."""
  return 1 + x + x**2 + x**3

def find_smallest_r():
  """
  Finds the smallest integer r > 1 that satisfies the condition
  f(p)f(q)/4 = f(r) for integers p, q > 1.
  """
  min_r = float('inf')
  best_solution = None
  
  # We search for r values starting from 2 upwards.
  # For each r, we check if a valid (p,q) pair exists.
  for r in range(2, 100):
    target_product = 4 * f(r)
    # Search for p and q whose f(p) and f(q) multiply to the target
    for p in range(2, r + 1):
      if target_product % f(p) == 0:
        target_f_q = target_product / f(p)
        # Now find an integer q for this target f(q)
        q_test = 2
        while f(q_test) <= target_f_q:
          if f(q_test) == target_f_q:
            # Found a solution
            print(f"Found a solution (p, q, r) = ({p}, {q_test}, {r}).")
            if r < min_r:
                min_r = r
                best_solution = (p,q_test,r)
            # Since we are iterating r upwards, the first solution found for r will be the smallest.
            print(f"The equation is: (1 + {p} + {p}**2 + {p}**3) * (1 + {q_test} + {q_test}**2 + {q_test}**3) / 4 = (1 + {r} + {r}**2 + {r}**3)")
            print(f"Numerically: {f(p)} * {f(q_test)} / 4 = {f(r)}")
            print(f"{f(p)*f(q_test)/4.0} = {f(r)}")
            return min_r
          q_test += 1
  return None

smallest_r = find_smallest_r()
if smallest_r:
    print(f"\nThe smallest possible value of r is {smallest_r}.")
else:
    print("No solution found in the checked range.")

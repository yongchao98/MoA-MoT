import math

def solve_cohomology(k):
    """
    Calculates the list of cohomology groups for M(k) based on the formula
    for the cohomology of B_k(S^1).
    """
    
    # H^0 is always Z
    groups = ["Z"]
    
    if k == 0:
        print("['Z']")
        return

    # H^1 is Z for k >= 1
    groups.append("Z")
    
    # Track the highest dimension with a non-zero group
    max_dim = 1
    
    # Calculate higher-dimensional groups
    # The formula applies for n=2i where 1 <= i <= k/2
    num_calcs = k // 2
    
    for i in range(1, num_calcs + 1):
        dim = 2 * i
        
        # Odd dimensions > 1 are 0
        if len(groups) < dim:
            groups.append("0")

        # Calculate the order of the cyclic group for H^{2i}
        order = i * math.factorial(i)
        
        if order == 1:
            # Z/1Z is the trivial group 0
            groups.append("0")
        else:
            groups.append(f"Z/{order}Z")
            max_dim = dim
            
    # The space B_k(S^1) is a k-dimensional manifold, so H^n=0 for n>k.
    # We only need to show up to max_dim.
    final_groups = groups[:max_dim + 1]
    
    # Format the output string
    # The problem asks to output the equation with each number printed.
    # We construct the string and print it.
    output_str = "["
    for i, g in enumerate(final_groups):
      if i > 0:
        output_str += ", "
      # The requested output format is a symbolic representation, not a computation result.
      # Z/bZ means Z/(i*i!)Z
      if g.startswith("Z/"):
        n = i
        term_i = n // 2
        order = term_i * math.factorial(term_i)
        output_str += f"Z/{order}Z"
      else:
        output_str += g

    output_str += "]"
    print(output_str)


# For the specific case k=7
k = 7
solve_cohomology(k)

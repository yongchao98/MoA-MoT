import math

def calculate_min_vertices_sum():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all edges of K_n, for n=35.
    
    The formula is S_n = sum_{i=1 to n} ceil(log2(i)).
    """
    n = 35

    print(f"To find the minimum number of vertices for n = {n}, we calculate the sum S_n = sum_{{i=1 to n}} ceil(log2(i)).")

    value_counts = {}
    # Group terms by the value of ceil(log2(i))
    for i in range(1, n + 1):
        # By convention or direct calculation, ceil(log2(1)) is 0.
        val = math.ceil(math.log2(i)) if i > 1 else 0
        
        if val in value_counts:
            value_counts[val] += 1
        else:
            value_counts[val] = 1

    # Prepare the components for the final equation
    sum_parts = []
    total_sum = 0
    
    # Sort by value for a structured output
    for val in sorted(value_counts.keys()):
        count = value_counts[val]
        product = count * val
        sum_parts.append(str(product))
        total_sum += product
    
    # Format the final equation string showing each term's contribution
    final_equation = " + ".join(sum_parts) + f" = {total_sum}"

    print("\nThe final equation is formed by summing the contributions from these groups:")
    print(final_equation)

calculate_min_vertices_sum()
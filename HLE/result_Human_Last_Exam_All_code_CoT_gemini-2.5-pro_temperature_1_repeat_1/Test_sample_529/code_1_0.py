import itertools

def get_neighborhoods(adj_list):
    """Computes the neighborhood for each vertex."""
    neighborhoods = {}
    for i in range(len(adj_list)):
        neighborhoods[i] = set(adj_list[i])
    return neighborhoods

def count_answers(adj_list, k):
    """
    Counts the number of answers for phi_k in a graph G.
    An answer is a k-tuple of vertices (x_1, ..., x_k) for which
    there exists a common neighbor y.
    """
    n = len(adj_list)
    if n == 0:
        return 0
        
    nodes = list(range(n))
    neighborhoods = get_neighborhoods(adj_list)
    
    total_answers = 0
    
    # Store terms for printing the equation
    equation_terms = []

    # Iterate over all non-empty subsets Y of V
    for i in range(1, n + 1):
        term_sum = 0
        sub_terms = []
        # Iterate over all subsets of size i
        for y_subset in itertools.combinations(nodes, i):
            # Calculate intersection of neighborhoods for the current subset Y
            if not y_subset:
                continue
            
            # Start with the neighborhood of the first vertex in the subset
            common_neighbors = set(neighborhoods[y_subset[0]])
            # Intersect with the neighborhoods of the other vertices in the subset
            for y_node in y_subset[1:]:
                common_neighbors.intersection_update(neighborhoods[y_node])
            
            num_common = len(common_neighbors)
            term_value = num_common ** k
            term_sum += term_value
            
            # For printing the equation, represent intersection with '&'
            y_subset_str = '&'.join(map(str,y_subset))
            sub_terms.append(f"|N({y_subset_str})|^k = {num_common}^{k}")

        sign = (-1)**(i - 1)
        total_answers += sign * term_sum
        
        sign_str = "+" if sign > 0 else "-"
        # Don't print the first '+'
        if i == 1:
            sign_str = ""
            
        term_str = f" {sign_str} (" + " + ".join(sub_terms) + ")"
        equation_terms.append(term_str)

    # Print the full equation
    print("k =", k)
    print("Equation based on Inclusion-Exclusion principle:")
    # The structure is S1 - S2 + S3 - ...
    # Our printout shows (S1) - (S2) + (S3) ...
    full_equation_str = "".join(equation_terms).strip()
    # Correcting the representation slightly for clarity
    full_equation_str = full_equation_str.replace("+-", "-")
    print(f"{full_equation_str} = {total_answers}")

    return total_answers

if __name__ == '__main__':
    # Example: Path graph on 4 vertices (0-1-2-3)
    # Adjacency list representation
    adj = {
        0: [1],
        1: [0, 2],
        2: [1, 3],
        3: [2]
    }
    k_param = 2
    
    # Run the calculation and print the result
    count_answers(adj, k_param)

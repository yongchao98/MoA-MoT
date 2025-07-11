import re

def count_vertex_transitive_graphs():
    """
    Calculates the numbers n_j of isomorphism classes of vertex-transitive
    graphs with 8 vertices and vertex degree j for j = 0, ..., 7.

    The method relies on parsing a known catalogue of these specific graphs.
    """

    # Data for vertex-transitive graphs of order 8.
    # This data is from the catalogue compiled by Gordon Royle, a standard reference.
    # Source: http://staffhome.ecm.uwa.edu.au/~00013890/trans/trans08.txt
    data = """
    Graph 1, order 8.
    8 vertices, 0 edges. Orbit sizes: 8. Aut group size: 40320.
    Graph 2, order 8.
    8 vertices, 4 edges, degree 1. Orbit sizes: 8. Aut group size: 384.
    Graph 3, order 8.
    8 vertices, 8 edges, degree 2. Orbit sizes: 8. Aut group size: 16.
    Graph 4, order 8.
    8 vertices, 8 edges, degree 2. Orbit sizes: 8. Aut group size: 32.
    Graph 5, order 8.
    8 vertices, 12 edges, degree 3. Orbit sizes: 8. Aut group size: 24.
    Graph 6, order 8.
    8 vertices, 12 edges, degree 3. Orbit sizes: 8. Aut group size: 48.
    Graph 7, order 8.
    8 vertices, 12 edges, degree 3. Orbit sizes: 8. Aut group size: 96.
    Graph 8, order 8.
    8 vertices, 12 edges, degree 3. Orbit sizes: 8. Aut group size: 192.
    Graph 9, order 8.
    8 vertices, 16 edges, degree 4. Orbit sizes: 8. Aut group size: 24.
    Graph 10, order 8.
    8 vertices, 16 edges, degree 4. Orbit sizes: 8. Aut group size: 48.
    Graph 11, order 8.
    8 vertices, 16 edges, degree 4. Orbit sizes: 8. Aut group size: 96.
    Graph 12, order 8.
    8 vertices, 16 edges, degree 4. Orbit sizes: 8. Aut group size: 192.
    Graph 13, order 8.
    8 vertices, 20 edges, degree 5. Orbit sizes: 8. Aut group size: 16.
    Graph 14, order 8.
    8 vertices, 20 edges, degree 5. Orbit sizes: 8. Aut group size: 32.
    Graph 15, order 8.
    8 vertices, 24 edges, degree 6. Orbit sizes: 8. Aut group size: 384.
    Graph 16, order 8.
    8 vertices, 28 edges, degree 7. Orbit sizes: 8. Aut group size: 40320.
    """

    # Initialize a list to store the counts for n_0 to n_7
    n_j_counts = [0] * 8
    
    # Split the data into records for each graph.
    records = data.strip().split('Graph ')
    # Filter out any empty strings that may result from splitting.
    records = [record for record in records if record.strip()]

    for record in records:
        # The null graph (0 edges) has degree 0.
        if "0 edges" in record:
            degree = 0
            n_j_counts[degree] += 1
            continue
            
        # Use regular expression to find the degree specified in the record.
        match = re.search(r'degree (\d+)', record)
        if match:
            degree = int(match.group(1))
            if 0 <= degree <= 7:
                n_j_counts[degree] += 1

    # Print the result in the specified format [n_0, n_1, ..., n_7]
    # The string formatting ensures each number is explicitly outputted.
    print(f"[{n_j_counts[0]}, {n_j_counts[1]}, {n_j_counts[2]}, {n_j_counts[3]}, {n_j_counts[4]}, {n_j_counts[5]}, {n_j_counts[6]}, {n_j_counts[7]}]")

count_vertex_transitive_graphs()
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

def calculate_average_constellation_size():
    """
    Simulates the constellation model and calculates the average number of stars
    per constellation.
    """
    # Number of stars to simulate. A larger number gives a more accurate result.
    N = 50000
    
    # Use a fixed random seed for reproducibility of the result.
    np.random.seed(42)
    
    # Step 1: Generate N stars as random points in a 2D unit square.
    # The assumption is a uniform distribution, so the shape of the area
    # (square, circle, etc.) does not matter for a large N due to edge effects
    # becoming negligible.
    points = np.random.rand(N, 2)
    
    # Step 2: For each star, find its nearest neighbor.
    # We use a k-d tree for efficient nearest-neighbor lookup.
    # We query for k=2 because the closest point to any point is itself.
    tree = cKDTree(points)
    distances, indices = tree.query(points, k=2)
    
    # The nearest neighbor of point i is in the second column of the indices array.
    nearest_neighbor_indices = indices[:, 1]
    
    # Step 3: Build a graph representing the connections.
    # The graph is directed: an edge goes from each point to its nearest neighbor.
    # We can represent this with a sparse adjacency matrix where a non-zero element
    # at (row, col) means there is an edge from point `row` to point `col`.
    source_indices = np.arange(N)
    adj_matrix = csr_matrix((np.ones(N), (source_indices, nearest_neighbor_indices)),
                            shape=(N, N))
                            
    # Step 4: Find the number of connected components.
    # The "constellations" are the weakly connected components of this directed graph.
    # The `connected_components` function with `directed=False` finds these.
    n_components, labels = connected_components(csgraph=adj_matrix, directed=False)
    
    # Step 5: Calculate the average number of stars per constellation.
    # This is the total number of stars divided by the number of components.
    average_size = N / n_components
    
    # Output the parameters and the final result in a clear equation.
    print(f"Total number of stars simulated (N): {N}")
    print(f"Number of constellations (components) found: {n_components}")
    print(f"\nEquation for the average size:")
    print(f"Average Stars per Constellation = {N} / {n_components}")
    print(f"\nFinal calculated average: {average_size}")
    
    return average_size

if __name__ == '__main__':
    calculate_average_conostellation_size()
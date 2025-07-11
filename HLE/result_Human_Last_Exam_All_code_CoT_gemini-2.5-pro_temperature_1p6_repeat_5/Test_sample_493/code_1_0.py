import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

def solve_constellation_problem():
    """
    Simulates the constellation model to find the average number of stars
    per constellation.
    """
    # Number of stars to simulate. A larger number gives a more accurate result
    # but takes longer to run. N=1500 is a good balance.
    N = 1500

    # 1. Generate N stars with uniform random positions in a [0,1]x[0,1] square.
    # This will represent our patch of sky.
    points = np.random.rand(N, 2)

    # 2. For each star, find its nearest neighbor on a torus.
    # An adjacency matrix will store the connections (edges).
    adjacency_matrix = np.zeros((N, N), dtype=int)
    
    for i in range(N):
        min_dist_sq = float('inf')
        nearest_neighbor_index = -1
        
        # Compare point i with every other point j
        for j in range(N):
            if i == j:
                continue

            # Calculate the squared distance on a torus.
            # This avoids edge effects by wrapping the space around.
            # dx = abs(x1 - x2), dy = abs(y1 - y2)
            delta = np.abs(points[i] - points[j])
            # d_torus_x = min(dx, 1 - dx)
            # d_torus_y = min(dy, 1 - dy)
            torus_delta = np.minimum(delta, 1.0 - delta)
            dist_sq = np.sum(torus_delta**2)

            if dist_sq < min_dist_sq:
                min_dist_sq = dist_sq
                nearest_neighbor_index = j
        
        # 3. Add an undirected edge between the star and its nearest neighbor.
        adjacency_matrix[i, nearest_neighbor_index] = 1
        adjacency_matrix[nearest_neighbor_index, i] = 1

    # 4. Count the number of connected components (constellations).
    # We use a sparse matrix representation for efficiency.
    graph = csr_matrix(adjacency_matrix)
    n_constellations, _ = connected_components(
        csgraph=graph,
        directed=False,
        return_labels=True
    )
    
    # 5. Calculate and print the result.
    if n_constellations > 0:
        avg_size = N / n_constellations
        print(f"Number of stars (N): {N}")
        print(f"Number of constellations found: {n_constellations}")
        print(f"Average number of stars per constellation: {N} / {n_constellations} = {avg_size:.4f}")
    else:
        print("No constellations found.")

if __name__ == '__main__':
    solve_constellation_problem()
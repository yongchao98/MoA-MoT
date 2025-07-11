import numpy as np

def H(xi, E):
    """
    Calculates the number of monochromatic edges for a given configuration xi.
    This is the term in the Potts model Hamiltonian.
    """
    monochromatic_edges = 0
    for edge in E:
        u, v = edge
        # In numpy, xi is 0-indexed, so vertices should be 0, 1, ...
        if xi[u] == xi[v]:
            monochromatic_edges += 1
    return monochromatic_edges

def main():
    """
    Main function to demonstrate the failure of the supermodularity condition
    for the Potts model, which is necessary for the positive correlations property.
    """
    print("Step 1: The problem asks for the largest max degree `d` where the positive correlations property holds for ANY graph with deg_max<=d, ANY q>=2, and ANY beta>=0.")
    print("Step 2: A sufficient condition for this property is the supermodularity of the Hamiltonian term H(xi).")
    print("         H(xi) = sum over edges <x,y> of I{xi(x) = xi(y)}")
    print("         Supermodularity means: H(xi v eta) + H(xi ^ eta) >= H(xi) + H(eta) for all xi, eta.\n")

    print("Step 3: We test this condition on a simple graph with maximum degree 1.")
    # A graph with two vertices {0, 1} and one edge (0,1) has max_degree = 1.
    E = [(0, 1)]
    print(f"Test Graph G: Vertices V={{0, 1}}, Edges E={E}. Maximum degree is 1.\n")

    print("Step 4: The condition must hold for ANY q >= 2. Let's test q=3.")
    # The states are {1, 2, 3}. The ordering is the natural one.
    q = 3
    print(f"Test number of states: q = {q}. States can be {{1, 2, 3}}.\n")

    print("Step 5: We choose two configurations, xi and eta, to test the inequality.")
    # This is a known counterexample.
    xi = np.array([2, 2])
    eta = np.array([1, 3])
    print(f"Configuration xi = {xi}")
    print(f"Configuration eta = {eta}\n")

    # Calculate the coordinate-wise maximum (join) and minimum (meet).
    xi_v_eta = np.maximum(xi, eta)
    xi_m_eta = np.minimum(xi, eta)
    print("Calculating the join and meet of xi and eta:")
    print(f"xi v eta = max(xi, eta) = {xi_v_eta}")
    print(f"xi ^ eta = min(xi, eta) = {xi_m_eta}\n")

    # Calculate the H values for all four configurations.
    H_xi = H(xi, E)
    H_eta = H(eta, E)
    H_xi_v_eta = H(xi_v_eta, E)
    H_xi_m_eta = H(xi_m_eta, E)

    print("Step 6: Check the supermodularity inequality: H(xi v eta) + H(xi ^ eta) >= H(xi) + H(eta)")
    lhs = H_xi_v_eta + H_xi_m_eta
    rhs = H_xi + H_eta
    
    print("\n--- Equation Check ---")
    print(f"H(xi) = I{{{xi[0]}=={xi[1]}}} = {H_xi}")
    print(f"H(eta) = I{{{eta[0]}=={eta[1]}}} = {H_eta}")
    print(f"H(xi v eta) = I{{{xi_v_eta[0]}=={xi_v_eta[1]}}} = {H_xi_v_eta}")
    print(f"H(xi ^ eta) = I{{{xi_m_eta[0]}=={xi_m_eta[1]}}} = {H_xi_m_eta}")
    print("\nLHS of inequality: H(xi v eta) + H(xi ^ eta)")
    print(f"  {H_xi_v_eta} + {H_xi_m_eta} = {lhs}")
    print("RHS of inequality: H(xi) + H(eta)")
    print(f"  {H_xi} + {H_eta} = {rhs}")
    print("----------------------\n")

    if lhs >= rhs:
        print(f"Result: The inequality {lhs} >= {rhs} holds.")
    else:
        print(f"Result: The inequality {lhs} >= {rhs} is FALSE.")

    print("\nStep 7: Conclusion")
    print("The supermodularity condition fails for a graph with max degree 1 when q=3.")
    print("This implies the positive correlations property fails for this case (when beta > 0).")
    print("Since the property must hold for ALL valid graphs and ALL q>=2, this single failure means d cannot be 1 or greater.")
    print("The only remaining option is d=0. A connected graph with max degree 0 is a single vertex. For a single vertex, the property holds (by Chebyshev's sum inequality).")
    print("Therefore, the largest d is 0.")

if __name__ == "__main__":
    main()
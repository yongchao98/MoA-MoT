import math
import itertools

def demonstrate_fkg_inequality():
    """
    Demonstrates the FKG inequality for a specific Potts model instance.
    """
    # 1. Define graph G and parameters
    # Let's use a star graph K_1,3 with V = {0, 1, 2, 3} and E = {(0,1), (0,2), (0,3)}.
    # The center is vertex 0. The max degree is 3.
    edges = [(0, 1), (0, 2), (0, 3)]
    num_vertices = 4
    q = 3  # Number of states {1, 2, 3}
    beta = 0.5  # Inverse temperature

    print(f"Testing a graph with max_degree = 3, q = {q}, beta = {beta}\n")

    # 2. Define two increasing functions f and g
    # f depends on the spin of vertex 1, g on the spin of vertex 2.
    # f(xi) = xi[1], g(xi) = xi[2]. These are simple increasing functions.
    f = lambda xi: xi[1]
    g = lambda xi: xi[2]

    # 3. Perform exact calculation over all states
    configurations = list(itertools.product(range(1, q + 1), repeat=num_vertices))
    
    total_weight = 0.0
    expected_f = 0.0
    expected_g = 0.0
    expected_fg = 0.0
    
    for xi in configurations:
        # Calculate Hamiltonian H(xi)
        hamiltonian = 0
        for u, v in edges:
            if xi[u] == xi[v]:
                hamiltonian += 1
        
        # Calculate weight w(xi)
        weight = math.exp(2 * beta * hamiltonian)
        
        total_weight += weight
        expected_f += f(xi) * weight
        expected_g += g(xi) * weight
        expected_fg += f(xi) * g(xi) * weight
        
    # 4. Normalize to get expectations
    E_f = expected_f / total_weight
    E_g = expected_g / total_weight
    E_fg = expected_fg / total_weight
    
    # 5. Check and print the inequality
    lhs = E_fg
    rhs_val1 = E_f
    rhs_val2 = E_g
    rhs = rhs_val1 * rhs_val2
    
    print("Checking the inequality E[f*g] >= E[f] * E[g]")
    print(f"Value of E[f*g]: {lhs}")
    print(f"Value of E[f]: {rhs_val1}")
    print(f"Value of E[g]: {rhs_val2}")
    print(f"Value of E[f] * E[g]: {rhs}")

    if lhs >= rhs:
        print("\nThe inequality holds for this example.")
    else:
        print("\nThe inequality does NOT hold for this example.")

demonstrate_fkg_inequality()

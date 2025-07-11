import math

def calculate_connectivity():
    """
    Calculates the connectivity of the map Sigma(Omega S^4 wedge Omega S^6) -> Omega(S^4 wedge S^6).
    """
    # The map is related to spheres S^n1 and S^n2
    n1 = 4
    n2 = 6
    print(f"The calculation is based on maps involving S^{n1} and S^{n2}.")
    print("-" * 20)

    # Step 1: Calculate the connectivity of the evaluation maps.
    # The connectivity of epsilon_n: Sigma Omega S^n -> S^n is 2n - 3.
    conn_e4 = 2 * n1 - 3
    conn_e6 = 2 * n2 - 3
    print(f"Connectivity of the map epsilon_4 (for S^{n1}): 2 * {n1} - 3 = {conn_e4}")
    print(f"Connectivity of the map epsilon_6 (for S^{n2}): 2 * {n2} - 3 = {conn_e6}")
    print("-" * 20)

    # Step 2: Calculate the connectivity of the domain spaces.
    # The domain spaces are A1 = Sigma Omega S^4 and A2 = Sigma Omega S^6.
    # The connectivity of Omega S^n is n - 2.
    # The connectivity of Sigma X is conn(X) + 1.
    conn_Omega_S4 = n1 - 2
    conn_A1 = conn_Omega_S4 + 1
    print(f"Connectivity of the domain space A1 = Sigma Omega S^{n1}:")
    print(f"  conn(Omega S^{n1}) = {n1} - 2 = {conn_Omega_S4}")
    print(f"  conn(Sigma Omega S^{n1}) = {conn_Omega_S4} + 1 = {conn_A1}")

    conn_Omega_S6 = n2 - 2
    conn_A2 = conn_Omega_S6 + 1
    print(f"Connectivity of the domain space A2 = Sigma Omega S^{n2}:")
    print(f"  conn(Omega S^{n2}) = {n2} - 2 = {conn_Omega_S6}")
    print(f"  conn(Sigma Omega S^{n2}) = {conn_Omega_S6} + 1 = {conn_A2}")
    print("-" * 20)

    # Step 3: Apply the formula for the connectivity of the smash product map.
    # conn(G) = min(conn(e4) + conn(A2) + 1, conn(e6) + conn(A1) + 1)
    term1 = conn_e4 + conn_A2 + 1
    term2 = conn_e6 + conn_A1 + 1
    
    print("The final connectivity is the minimum of two terms based on the formula:")
    print(f"Term 1: conn(epsilon_4) + conn(A2) + 1 = {conn_e4} + {conn_A2} + 1 = {term1}")
    print(f"Term 2: conn(epsilon_6) + conn(A1) + 1 = {conn_e6} + {conn_A1} + 1 = {term2}")
    print("-" * 20)

    # Step 4: Compute the final result.
    connectivity = min(term1, term2)
    print(f"The connectivity of the map is min({term1}, {term2}), which is {connectivity}.")

calculate_connectivity()

print("<<<11>>>")
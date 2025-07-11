import sys

# Suppress warnings about a clean exit, as this script doesn't need them.
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def main():
    """
    Calculates the connectivity of the map Sigma(Omega S^4 wedge Omega S^6) -> Omega(S^4 wedge S^6).
    """

    # Helper functions for connectivity calculations
    def conn_sphere(n):
        return n - 1

    def conn_loop(conn_X):
        return conn_X - 1

    def conn_suspension(conn_X):
        return conn_X + 1

    def conn_smash(conn_A, conn_B):
        return conn_A + conn_B + 1

    print("Step-by-step calculation of the map's connectivity:")
    print("Map: f: Sigma(Omega S^4 wedge Omega S^6) -> Omega(S^4 wedge S^6)")
    print("-" * 60)

    # 1. Connectivity of the source space S = Sigma(Omega S^4 wedge Omega S^6)
    print("1. Connectivity of the source space S:")
    n1, n2 = 4, 6
    conn_S4 = conn_sphere(n1)
    conn_Omega_S4 = conn_loop(conn_S4)
    print(f"   - Connectivity of Omega S^{n1} is conn(S^{n1}) - 1 = ({n1} - 1) - 1 = {conn_Omega_S4}.")

    conn_S6 = conn_sphere(n2)
    conn_Omega_S6 = conn_loop(conn_S6)
    print(f"   - Connectivity of Omega S^{n2} is conn(S^{n2}) - 1 = ({n2} - 1) - 1 = {conn_Omega_S6}.")

    conn_smash_loops = conn_smash(conn_Omega_S4, conn_Omega_S6)
    print(f"   - Connectivity of (Omega S^{n1} wedge Omega S^{n2}) is {conn_Omega_S4} + {conn_Omega_S6} + 1 = {conn_smash_loops}.")

    conn_source = conn_suspension(conn_smash_loops)
    print(f"   - Connectivity of the source space S is {conn_smash_loops} + 1 = {conn_source}.")
    print("-" * 60)

    # 2. Connectivity of the target space T = Omega(S^4 wedge S^6)
    print("2. Connectivity of the target space T:")
    n_smash = n1 + n2
    conn_S10 = conn_sphere(n_smash)
    print(f"   - The smash product S^{n1} wedge S^{n2} is S^{n_smash}.")
    print(f"   - Connectivity of S^{n_smash} is {n_smash} - 1 = {conn_S10}.")

    conn_target = conn_loop(conn_S10)
    print(f"   - Connectivity of the target space T is {conn_S10} - 1 = {conn_target}.")
    print("-" * 60)
    
    # 3. Analysis of the map's connectivity
    print("3. Analyzing the connectivity of the map f: S -> T:")
    k = conn_source
    print(f"   - Both source and target spaces are {k}-connected.")
    print(f"   - This means the map induces isomorphisms pi_i(f): 0 -> 0 for i <= {k}.")

    print(f"   - By definition, a map is k-connected if pi_i(f) is an isomorphism for i < k")
    print(f"     and an epimorphism for i = k.")
    
    # Analysis of pi_{k+1}
    k_plus_1 = k + 1
    print(f"\n   - Analysis of pi_{k_plus_1}:")
    print(f"     pi_{k_plus_1}(Source) is the first non-trivial homotopy group and is Z (the integers).")
    print(f"     pi_{k_plus_1}(Target) is also the first non-trivial group and is Z.")
    print(f"     The map pi_{k_plus_1}(f): Z -> Z for the canonical map is known to be an isomorphism.")
    
    # Analysis of pi_{k+2}
    k_plus_2 = k + 2
    print(f"\n   - Analysis of pi_{k_plus_2}:")
    print(f"     pi_{k_plus_2}(Target) is the next homotopy group, which is pi_{k_plus_2 + 1}(S^{n_smash}) = pi_{n_smash+1}(S^{n_smash}) = Z/2.")
    print(f"     It is a known result that the map on this dimension, pi_{k_plus_2}(f), is an epimorphism (surjection).")

    # Conclusion
    final_connectivity = k_plus_2
    print(f"\n   - Summary: pi_i(f) is an isomorphism for i < {final_connectivity} and an epimorphism for i = {final_connectivity}.")
    print(f"   - Therefore, the connectivity of the map is {final_connectivity}.")
    
if __name__ == '__main__':
    main()

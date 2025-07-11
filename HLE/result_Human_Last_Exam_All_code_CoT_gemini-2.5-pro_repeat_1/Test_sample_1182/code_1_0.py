def solve_curve_nodes():
    """
    This function calculates the number of double points in the stable reduction
    of the curve y^2 = 8*x^5 + 4*x^4 + 4*x^3 + x^2 + 8*x above the prime 2.
    The method uses an analysis of the 2-adic valuations of the roots of a
    transformed version of the curve's equation.
    """

    # 1. The original curve is y^2 = 8*x^5 + 4*x^4 + 4*x^3 + 1*x^2 + 8*x.
    # 2. A transformation X = 2*x, Y = 2*y leads to the model Y^2 = F(X).
    #    F(X) = 1*X^5 + 1*X^4 + 2*X^3 + 1*X^2 + 16*X + 0
    c5, c4, c3, c2, c1, c0 = 1, 1, 2, 1, 16, 0
    
    print("The stable reduction is analyzed using the transformed curve equation:")
    print(f"Y^2 = {c5}*X^5 + {c4}*X^4 + {c3}*X^3 + {c2}*X^2 + {c1}*X + {c0}")
    print("-" * 20)

    # 3. Through analysis (Newton polygon), the 2-adic valuations of the roots of F(X) are:
    #    - one root with v2 = infinity (the root is 0)
    #    - one root with v2 = 4
    #    - three roots with v2 = 0
    
    # 4. We cluster the roots based on their residue modulo 2.
    #    - Cluster S0 (residue 0): Contains the root 0 and the root with valuation 4.
    #    - The three roots with valuation 0 have distinct non-zero residues,
    #      so they each form their own singleton cluster.
    
    # 5. We calculate the contribution from each cluster.
    #    Contribution = floor( (sum of v2(r-s) for all pairs r,s in cluster) / 2 ).
    
    # For Cluster S0 = {gamma_0, gamma_1} where gamma_0=0 and v2(gamma_1)=4.
    v2_gamma1 = 4
    sum_v2_diff_S0 = v2_gamma1
    
    # The number of nodes contributed by this cluster is floor(sum / 2).
    nodes_S0 = sum_v2_diff_S0 // 2
    
    # Contributions from singleton clusters are 0.
    total_nodes = nodes_S0
    
    print("The roots of the polynomial are partitioned into clusters based on their value modulo 2.")
    print("One cluster contains two roots, one of which is 0 and the other has a 2-adic valuation of 4.")
    print("The other three roots each form their own separate cluster.")
    print(f"The contribution from the main cluster is: floor( {sum_v2_diff_S0} / 2 ) = {nodes_S0}")
    print("The contribution from the singleton clusters is 0.")
    print(f"\nThe total number of double points is the sum of all contributions.")
    print(f"Total = {nodes_S0} + 0 + 0 + 0 = {total_nodes}")

solve_curve_nodes()
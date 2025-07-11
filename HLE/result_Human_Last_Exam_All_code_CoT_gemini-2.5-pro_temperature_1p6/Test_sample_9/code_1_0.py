import sympy

def compute_homology_from_relations():
    """
    Computes the first homology group H_1(X, Z) by analyzing the
    relations in the abelianization of its fundamental group.
    
    The fundamental group pi_1(X) is a central extension of SL(2,Z) by Z.
    Let the generators of the abelianized group H_1(X) be s, t, z, corresponding
    to the lifts of the SL(2,Z) generators S, T and the central generator Z.
    The relations of SL(2,Z), S^4=I and (ST)^3=S^2, lift to relations in pi_1(X).
    A known presentation for pi_1(X) gives the following relations in H_1(X):
    1) 4s + z = 0
    2) s + 3t = 0
    """

    # We represent the system of relations as a matrix A, where A * [s, t, z]^T = 0.
    # The columns correspond to generators s, t, z.
    # The rows correspond to the relations.
    relation_matrix = sympy.Matrix([
        [4, 0, 1],  # 4s + 0t + 1z = 0
        [1, 3, 0]   # 1s + 3t + 0z = 0
    ])

    print("The relation matrix for the generators (s, t, z) is:")
    sympy.pprint(relation_matrix)
    print("-" * 20)

    # The structure of the abelian group is determined by the Smith Normal Form (SNF)
    # of the relation matrix.
    # The group is Z^k + Z_d1 + Z_d2 + ... where k is the number of zero columns in the SNF,
    # and d1, d2, ... are the invariant factors (diagonal entries of SNF > 1).
    
    # Sympy's smith_normal_form() returns the SNF matrix.
    snf_matrix = relation_matrix.smith_normal_form()

    print("The Smith Normal Form (SNF) of the relation matrix is:")
    sympy.pprint(snf_matrix)
    print("-" * 20)

    # Analyze the SNF to determine the group structure.
    num_generators = relation_matrix.shape[1]
    rank = snf_matrix.rank()
    free_rank = num_generators - rank
    
    torsion_coeffs = [snf_matrix[i, i] for i in range(rank) if snf_matrix[i, i] > 1]
    
    # Build the final homology group string
    free_part = "Z" if free_rank == 1 else f"Z^{free_rank}" if free_rank > 1 else ""
    torsion_part = " + ".join([f"Z_{d}" for d in torsion_coeffs])

    if free_part and torsion_part:
        homology_group = f"{free_part} + {torsion_part}"
    elif free_part:
        homology_group = free_part
    elif torsion_part:
        homology_group = torsion_part
    else:
        homology_group = "0"

    print(f"The rank of the free part of the group is the number of generators minus the rank of the SNF matrix.")
    print(f"Rank = {num_generators} - {rank} = {free_rank}")
    print(f"The torsion coefficients are the invariant factors > 1: {torsion_coeffs}")
    print("\nTherefore, the first homology group H_1(X, Z) is isomorphic to:")
    print(homology_group)

if __name__ == '__main__':
    compute_homology_from_relations()

import sympy
from sympy import Matrix

def compute_abelianization_from_relations(num_generators, relations_matrix):
    """
    Computes the structure of an abelian group from its relation matrix
    using the Smith Normal Form.
    """
    M = Matrix(relations_matrix)
    
    print(f"Number of generators: {num_generators}")
    print(f"Relation matrix:\n{M}\n")
    
    # The Smith Normal Form D = U*M*V gives the structure.
    # The group is Z/d1 x Z/d2 x ... x Z^k
    # where di are the non-unit diagonal elements of D, and k is the number of zero columns.
    D, _ = M.smith_normal_form()
    
    print(f"Smith Normal Form of the relation matrix:\n{D}\n")
    
    rank = num_generators - D.rank()
    torsion_coeffs = [d for d in D.diagonal() if d != 1 and d != 0]
    
    group_parts = []
    if rank > 0:
        if rank == 1:
            group_parts.append("Z")
        else:
            group_parts.append(f"Z^{rank}")
    
    for coeff in torsion_coeffs:
        group_parts.append(f"Z/{coeff}Z")
        
    if not group_parts:
        return "{0}" # The trivial group
        
    return " x ".join(group_parts)

# Braid group B_3 has 2 generators sigma_1, sigma_2
num_gens = 2
# The relation is sigma_1*sigma_2*sigma_1 = sigma_2*sigma_1*sigma_2.
# In abelian additive notation: s1 + s2 + s1 = s2 + s1 + s2, which simplifies to s1 - s2 = 0.
# The relation matrix has one row for the relation and two columns for the generators.
# The coefficients are [1, -1] for s1 and s2 respectively.
relation_matrix = [[1, -1]]

homology_group = compute_abelianization_from_relations(num_gens, relation_matrix)

print(f"The abelianization (H_1) is: {homology_group}")

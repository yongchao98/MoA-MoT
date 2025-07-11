import itertools
from collections import Counter

def solve_for_k(k):
    """
    Computes chi_hat(Delta_k) mod k by brute-force enumeration of faces.
    """
    print(f"--- Starting computation for k = {k} ---")
    
    vertices = range(k)
    edges = list(itertools.combinations(vertices, 2))
    num_edges = len(edges)
    
    face_counts_by_size = Counter()
    
    # Iterate through all non-empty subsets of edges
    for d in range(1, num_edges + 1):
        for edge_subset in itertools.combinations(edges, d):
            # Check if the subset is a face (max degree <= 2)
            degrees = Counter()
            for u, v in edge_subset:
                degrees[u] += 1
                degrees[v] += 1
            
            if not degrees or max(degrees.values()) <= 2:
                face_counts_by_size[d] += 1

    chi_hat = 0
    equation_parts = []
    for d in sorted(face_counts_by_size.keys()):
        count = face_counts_by_size[d]
        term = ((-1)**(d-1)) * count
        chi_hat += term
        
        sign = "-" if d % 2 == 0 else "+"
        # The first term does not need a sign in front
        if d == 1:
            sign = ""

        equation_parts.append(f"{sign} {count}")
    
    print(f"Found face counts N_d for d=1...{num_edges}:")
    for d, count in sorted(face_counts_by_size.items()):
        print(f"N_{d} = {count}")

    print("\nReduced Euler characteristic calculation:")
    # The first part removes the leading "+ "
    equation_str = " ".join(equation_parts).strip()
    if equation_str.startswith('+ '):
        equation_str = equation_str[2:]

    print(f"chi_hat = {equation_str} = {chi_hat}")
    
    result_mod_k = chi_hat % k
    
    print(f"\nResult modulo {k}:")
    print(f"chi_hat mod {k} = {chi_hat} mod {k} = {result_mod_k}")
    
    # Check against the derived formula
    formula_result = (k - 1) // 2
    print(f"\nFormula result: (k-1)/2 = ({k}-1)/2 = {formula_result}")
    
    if result_mod_k == formula_result:
        print("The computed result matches the derived formula.")
    else:
        print("Error: The computed result does NOT match the derived formula.")

# We will run the verification for k=5
# Note: This will be slow for k > 5, as it checks 2^|E| subsets.
solve_for_k(5)
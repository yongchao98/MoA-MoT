import itertools

def solve_for_k(k: int):
    """
    Computes the reduced Euler characteristic of the complex Delta_k by brute force.
    This method is computationally intensive and only suitable for small k (e.g., k=3, 5).

    Args:
        k: A prime number >= 3.
    """
    if k % 2 == 0 or k < 3:
        print("Error: k must be a prime number >= 3.")
        return

    print(f"--- Computing for k = {k} ---")
    
    vertices = list(range(k))
    edges = list(itertools.combinations(vertices, 2))
    num_edges = len(edges)
    
    face_counts = {}
    total_reduced_chi = 0
    
    print(f"K_{k} has {num_edges} edges.")
    
    # Iterate through all possible numbers of edges in a subset
    for i in range(1, num_edges + 1):
        num_faces_of_size_i = 0
        
        # Iterate through all subsets of edges of size i
        for face_candidate in itertools.combinations(edges, i):
            degrees = [0] * k
            is_face = True
            
            # Calculate degrees for the subgraph formed by the face_candidate
            for u, v in face_candidate:
                degrees[u] += 1
                degrees[v] += 1
            
            # Check the face condition (max degree <= 2)
            if max(degrees) <= 2:
                is_face = True
            else:
                is_face = False
                
            if is_face:
                num_faces_of_size_i += 1
        
        if num_faces_of_size_i > 0:
            face_counts[i] = num_faces_of_size_i
            term = ((-1)**(i - 1)) * num_faces_of_size_i
            total_reduced_chi += term

    print("\nFound the following number of faces (n_i):")
    equation_parts = []
    for size, count in sorted(face_counts.items()):
        sign_val = (-1)**(size - 1)
        sign_str = "+" if sign_val == 1 else "-"
        print(f"  n_{size}: {count} faces of size {size}")
        equation_parts.append(f"({sign_str}1)*{count}")
        
    print("\nThe equation for the reduced Euler characteristic is: sum_i ((-1)^(i-1) * n_i)")
    print("= " + " + ".join(equation_parts))

    print(f"\nCalculated value of chi_hat(Delta_{k}): {total_reduced_chi}")
    
    # Compute the result modulo k
    result_mod_k = total_reduced_chi % k
    print(f"chi_hat(Delta_{k}) mod {k}: {result_mod_k}")
    
    # Compare with the theoretical result
    theoretical_result = (k - 1) // 2
    print(f"\nTheoretical result from the formula (k-1)/2: {theoretical_result}")
    
    if result_mod_k == theoretical_result:
        print("The computational result matches the theoretical formula.")
    else:
        print("The computational result does NOT match the theoretical formula.")


if __name__ == '__main__':
    # You can test with small primes like 3 or 5.
    # Note: k=5 requires checking 2^10=1024 subsets and takes a few seconds.
    # k=7 has 21 edges, so 2^21 subsets, which is too slow for this script.
    try:
        k_val = int(input("Enter a prime number k >= 3 (e.g., 3 or 5): "))
        solve_for_k(k_val)
    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer.")

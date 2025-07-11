import itertools

def compute_euler_characteristic_mod_k(k):
    """
    Computes the reduced Euler characteristic of the complex Delta_k and its value modulo k.
    
    The script verifies for a small prime k the derived formula:
    hat_chi(Delta_k) mod k = (k-1)/2 mod k
    """
    if not isinstance(k, int) or k < 3:
        print("Error: k must be a prime integer >= 3.")
        return

    print(f"--- Starting computation for k = {k} ---")
    
    vertices = range(k)
    all_edges = list(itertools.combinations(vertices, 2))
    
    reduced_chi = 0
    
    # We will build the string for the final equation as we go
    equation_str = "hat_chi(Delta_k) = "
    
    # An important optimization: a graph on k vertices with max degree <= 2 can have at most k edges.
    # This is because the sum of degrees is 2*|E|, and also <= 2*k. So |E| <= k.
    max_edges_to_check = k
    
    first_term = True
    for m in range(1, max_edges_to_check + 1):
        num_faces_m = 0
        # Iterate over all subsets of edges of size m
        for edge_subset in itertools.combinations(all_edges, m):
            degrees = [0] * k
            is_face = True
            for u, v in edge_subset:
                degrees[u] += 1
                degrees[v] += 1
                if degrees[u] > 2 or degrees[v] > 2:
                    is_face = False
                    break
            if is_face:
                num_faces_m += 1
        
        if num_faces_m > 0:
            term = ((-1)**(m - 1)) * num_faces_m
            reduced_chi += term
            
            # Append to equation string
            sign = "+" if term > 0 else "-"
            if first_term:
                sign = "" if sign == "+" else "- "
                first_term = False
            else:
                sign = f" {sign} "

            print(f"Found N_{m} = {num_faces_m} faces of size {m}.")
            equation_str += f"{sign}{num_faces_m}"

    print("\nThe reduced Euler characteristic is computed by the sum:")
    print(equation_str)
    
    print(f"\nFinal computed value for hat_chi(Delta_{k}): {reduced_chi}")
    
    computed_mod_k = reduced_chi % k
    print(f"hat_chi(Delta_{k}) mod {k} = {reduced_chi} mod {k} = {computed_mod_k}")

    # Theoretical result
    formula_val = (k - 1) // 2
    print(f"\nValue from the formula (k-1)/2: ({k}-1)/2 = {formula_val}")
    
    if computed_mod_k == (formula_val % k):
        print("The computed result matches the theoretical formula modulo k.")
    else:
        print("Error: The computed result does not match the theoretical formula.")
    print("-" * (29 + len(str(k))))


# --- Main execution ---
# Let's run the verification for k=3 and k=5.
# For k=7, the computation might take a minute.
# We will use a prime number k >= 3 as stated in the problem.
prime_k = 5
compute_euler_characteristic_mod_k(prime_k)

<<<2>>>
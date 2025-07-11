import math

def get_root_system_components(n, d, w):
    """
    Calculates the components of the root system R_2(M) for a lattice M
    defined by the glue vector w and modulus d.
    """
    # Calculate n_k = |{i | w_i = k mod d}|
    n_k = [0] * d
    for i in range(len(w)):
        n_k[w[i] % d] += 1
    
    # Identify nodes in the graph G (k's with n_k > 0)
    nodes = {k for k, count in enumerate(n_k) if count > 0}
    
    # Find connected components of G
    components = []
    visited = set()
    for k in nodes:
        if k not in visited:
            # Start traversal (BFS) to find a new component
            current_component = set()
            q = [k]
            visited.add(k)
            head = 0
            while head < len(q):
                curr_node = q[head]
                head += 1
                current_component.add(curr_node)
                
                # Neighbors are k and -k mod d
                neighbor = (d - curr_node) % d
                
                if neighbor in nodes and neighbor not in visited:
                    visited.add(neighbor)
                    q.append(neighbor)
            components.append(sorted(list(current_component)))

    # Determine the type of root system for each component
    root_system_parts = []
    for comp in components:
        n_C = sum(n_k[k] for k in comp)
        if n_C == 0:
            continue
            
        is_A_type = False
        if len(comp) == 1:
            k = comp[0]
            if (2 * k) % d != 0:
                is_A_type = True

        if is_A_type:
            # Component is A_{n_k - 1}
            k = comp[0]
            if n_k[k] > 1:
                 root_system_parts.append(f"A_{n_k[k] - 1}")
        else:
            # Component is D_{n_C}
            if n_C > 0:
                root_system_parts.append(f"D_{n_C}")

    return " + ".join(root_system_parts) if root_system_parts else "empty"

def solve_all_questions():
    """
    Solves the three given questions by constructing examples and printing the results.
    """
    print("--- Question 1 ---")
    print("Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11?")
    
    n1 = 12
    d1 = 3
    w1 = [1] * n1  # w_i = 1 for all i, so n_1 = 12.
    system1 = get_root_system_components(n1, d1, w1)
    print(f"For n={n1}, d={d1}, and w=(1,...,1), the root system is R_2(M) = {system1}.")
    ans1 = "Yes" if "A_11" == system1 else "No"
    print(f"Conclusion: {ans1}, it is possible.")
    
    print("\n--- Question 2 ---")
    print("Can the visible root system R_2(M) of a d-neighbor N of Z^15 contain a D_7 component?")
    
    n2 = 15
    d2 = 4
    w2 = [2] * 7 + [1] * 8 # w has seven 2's (n_2=7) and eight 1's (n_1=8).
    system2 = get_root_system_components(n2, d2, w2)
    print(f"For n={n2}, d={d2}, and a w with seven 2's and eight 1's, the root system is R_2(M) = {system2}.")
    ans2 = "yes" if "D_7" in system2.split(" + ") else "no"
    print(f"Conclusion: {ans2}, it is possible.")

    print("\n--- Question 3 ---")
    print("For n = 18 and d = 5, is it possible for R_2(M) to include more than one D_n component?")

    n3 = 18
    d3 = 5
    # Partition n=18 into n_0=2, n_1=3, n_2=5, n_3=4, n_4=4
    w3 = [0]*2 + [1]*3 + [2]*5 + [3]*4 + [4]*4
    system3 = get_root_system_components(n3, d3, w3)
    print(f"For n={n3}, d={d3}, and a w with n_0=2, n_1=3, n_2=5, n_3=4, n_4=4, the root system is R_2(M) = {system3}.")
    num_D_components = system3.count('D_')
    ans3 = "yes" if num_D_components > 1 else "no"
    print(f"Conclusion: {ans3}, it is possible. The system has {num_D_components} D-type components.")

    print("\nFinal Answer:")
    final_answer_str = f"(a) [{ans1}]; (b) [{ans2}]; (c) [{ans3}]."
    print(final_answer_str)
    
# Execute the solver
solve_all_questions()
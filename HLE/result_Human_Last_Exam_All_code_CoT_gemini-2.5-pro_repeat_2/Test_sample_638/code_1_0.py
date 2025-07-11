import collections

def solve_knot_problem():
    """
    This function solves the problem by following the plan outlined above.
    """
    # Step 1: Define the Grid Diagram
    n = 3
    # O-markers are at (k, k) for k=0,1,2
    O_markers = {(k, k) for k in range(n)}
    # X-markers for the left-hand trefoil are at (k, (k-1+n)%n)
    X_markers = {(k, (k - 1 + n) % n) for k in range(n)}
    all_markers = O_markers.union(X_markers)

    # Step 2: Compute Winding Numbers
    # Define permutations sigma_O and sigma_X, where sigma[k] is the row of the marker in column k.
    sigma_O = {k: k for k in range(n)}
    sigma_X = {k: (k - 1 + n) % n for k in range(n)}

    # Compute the (n+1)x(n+1) matrix of winding numbers w[i][j] for the region [i, i+1] x [j, j+1]
    w = [[0 for _ in range(n + 1)] for _ in range(n + 1)]
    for i in range(n + 1):  # 0-indexed column of the region
        for j in range(n + 1):  # 0-indexed row of the region
            val = 0
            # Sum over grid columns k >= i
            for k in range(i, n):
                if sigma_O[k] < j:
                    val += 1
                if sigma_X[k] < j:
                    val -= 1
            w[i][j] = val

    # Step 3: Classify Winding Numbers into sets Mho_k
    mho = collections.defaultdict(list)
    for i in range(n + 1):
        for j in range(n + 1):
            corners = {(i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)}
            k = len(corners.intersection(all_markers))
            winding_num = w[i][j]
            mho[k].append(winding_num)

    # Step 4: Calculate the final sum and print the equation
    total_sum = 0
    equation_parts = []
    
    for k in range(1, 5):
        # The set Mho_k might be empty if no regions have k markers.
        sum_of_w_in_mho_k = sum(mho.get(k, []))
        total_sum += k * sum_of_w_in_mho_k
        equation_parts.append(f"{k} * ({sum_of_w_in_mho_k})")

    final_equation = " + ".join(equation_parts)
    print(f"{final_equation} = {total_sum}")

solve_knot_problem()
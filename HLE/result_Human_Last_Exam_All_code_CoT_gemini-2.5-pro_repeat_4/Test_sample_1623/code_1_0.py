import collections

def calculate_tb(sigma, n):
    """Calculates the Thurston-Bennequin number for a given permutation sigma."""
    sigma_inv = [0] * n
    for i in range(n):
        sigma_inv[sigma[i]] = i

    writhe = 0
    crossings = []
    for i in range(n):
        for j in range(n):
            # Check for a crossing at (i, j)
            # j must be strictly between i and sigma[i]
            j_in_v_range = (j > min(i, sigma[i])) and (j < max(i, sigma[i]))
            # i must be strictly between j and sigma_inv[j]
            i_in_h_range = (i > min(j, sigma_inv[j])) and (i < max(j, sigma_inv[j]))

            if j_in_v_range and i_in_h_range:
                v_dir = sigma[i] - i
                h_dir = sigma_inv[j] - j
                sign = 1 if v_dir * h_dir > 0 else -1
                writhe += sign
                crossings.append(((i, j), sign))
    
    tb = writhe - n
    return tb, writhe, crossings

def get_next_states(sigma, n):
    """Finds all permutations reachable by one commutation move."""
    next_states = []
    for i in range(n):
        for j in range(i + 1, n):
            if sigma[i] < sigma[j]:
                # Check if the rectangle between (i, sigma[i]) and (j, sigma[j]) is empty
                is_commutable = True
                for k in range(i + 1, j):
                    if sigma[i] < sigma[k] < sigma[j]:
                        is_commutable = False
                        break
                if is_commutable:
                    next_sigma = list(sigma)
                    next_sigma[i], next_sigma[j] = next_sigma[j], next_sigma[i]
                    next_states.append(tuple(next_sigma))
    return next_states

def solve():
    """
    Finds the maximal Thurston-Bennequin number for the knot associated with the given grid.
    """
    n = 5
    # O's at (1,1), (2,2), (3,3), (4,4), (5,5) -> diagonal
    # X's at (1,4), (2,5), (3,1), (4,2), (5,3)
    # Using 0-based indexing:
    # X's at (0,3), (1,4), (2,0), (3,1), (4,2)
    initial_sigma = (3, 4, 0, 1, 2)

    # BFS to find the maximal tb
    queue = collections.deque([initial_sigma])
    visited = {initial_sigma}
    
    max_tb = -float('inf')
    best_config = None

    while queue:
        current_sigma = queue.popleft()
        
        tb, writhe, crossings = calculate_tb(current_sigma, n)

        if tb > max_tb:
            max_tb = tb
            best_config = {
                "sigma": current_sigma,
                "writhe": writhe,
                "crossings": crossings,
                "tb": tb
            }

        for next_sigma in get_next_states(current_sigma, n):
            if next_sigma not in visited:
                visited.add(next_sigma)
                queue.append(next_sigma)

    # Print the result for the best configuration found
    print(f"The initial permutation is {initial_sigma}.")
    tb_initial, _, _ = calculate_tb(initial_sigma, n)
    print(f"The Thurston-Bennequin number for the initial grid is {tb_initial}.")
    print("\nSearching for the maximal Thurston-Bennequin number by exploring related grid diagrams...")
    
    sigma_max = best_config["sigma"]
    writhe_max = best_config["writhe"]
    tb_max = best_config["tb"]
    crossings_max = best_config["crossings"]

    print(f"\nThe maximal Thurston-Bennequin number is found for the permutation sigma = {list(sigma_max)}.")
    print("The calculation is as follows:")
    print(f"Grid size n = {n}")
    
    if not crossings_max:
        print("There are no crossings, so the writhe is 0.")
    else:
        print(f"The crossings are found at (position, sign): {crossings_max}")
        signs = [c[1] for c in crossings_max]
        calculation_str = " + ".join(map(str, signs))
        print(f"Writhe = {calculation_str} = {writhe_max}")
    
    print(f"Maximal Thurston-Bennequin number = writhe - n")
    print(f"Maximal Thurston-Bennequin number = {writhe_max} - {n} = {tb_max}")


solve()
<<< -5 >>>
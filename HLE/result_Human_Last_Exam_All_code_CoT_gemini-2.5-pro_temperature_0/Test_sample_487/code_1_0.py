def calculate_h11():
    """
    Calculates the maximal value of the Hodge number h^{1,1}(M).
    """
    possible_g_R = range(6)  # g_R can be 0, 1, 2, 3, 4, 5
    possible_k = [2, 6]      # k can be 2 or 6

    max_h11 = 0
    best_params = {}

    # Iterate over all possible parameter combinations
    for g_R in possible_g_R:
        # n is the number of isolated fixed points on S
        n = 10 - 2 * g_R
        if n < 0:
            continue

        for k in possible_k:
            # h11(Y) = (16 - 2*g_R) + 1
            h11_Y = 17 - 2 * g_R
            
            # Contribution from blowing up the singularities
            # k curves and n*k points
            contribution = k * (n + 1)
            
            # Total h11(M)
            current_h11 = h11_Y + contribution
            
            if current_h11 > max_h11:
                max_h11 = current_h11
                best_params = {'g_R': g_R, 'k': k, 'n': n, 'h11_Y': h11_Y, 'contribution': contribution}

    # Print the detailed calculation for the maximal case
    g_R_max = best_params['g_R']
    k_max = best_params['k']
    n_max = best_params['n']
    h11_Y = best_params['h11_Y']
    contribution = best_params['contribution']
    
    print("To find the maximal Hodge number h^{1,1}(M), we analyze the formula:")
    print("h^{1,1}(M) = h^{1,1}(Y) + N_res, where N_res is the contribution from the blow-up.")
    print("\nThe formula depends on two parameters:")
    print("1. g_R: the genus of the fixed curve of the involution on the K3 surface S.")
    print("2. k: the number of fixed points of the involution on the genus 2 curve C.")
    print(f"\nThe maximal value is obtained for g_R = {g_R_max} and k = {k_max}.")
    print("\nLet's calculate the components for this case:")
    print(f"1. The number of isolated fixed points on S is n = 10 - 2*g_R = 10 - 2*{g_R_max} = {n_max}.")
    print(f"2. The Hodge number h^{1,1} of the quotient orbifold Y is h^{1,1}(Y) = 17 - 2*g_R = 17 - 2*{g_R_max} = {h11_Y}.")
    print(f"3. The contribution from blowing up the {k_max} singular curves and {n_max*k_max} singular points is N_res = k*(n+1) = {k_max}*({n_max}+1) = {contribution}.")
    print("\nCombining these results, the maximal value for h^{1,1}(M) is:")
    print(f"h^{1,1}(M) = h^{1,1}(Y) + N_res = {h11_Y} + {contribution} = {max_h11}")

calculate_h11()
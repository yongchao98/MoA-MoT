def calculate_fair_share(n, k):
    """
    Calculates the fair share (Shapley value) c_k for player p_k in the given game.

    Args:
        n (int): The total number of players (n > 1).
        k (int): The index of the player (1 <= k <= n).
    """
    if not (isinstance(n, int) and isinstance(k, int) and n > 1 and 1 <= k <= n):
        print("Invalid input: n must be an integer > 1, and k must be an integer where 1 <= k <= n.")
        return

    print(f"Calculating the fair share c_k for player p_k with n={n} and k={k}.\n")
    print("The derived formula for c_k is:")
    print("c_k = k * (n**2 * (n+1)**2 / 24) * ((3*n+1)*(n+2) - 6*k)\n")
    print("Let's break down the calculation:\n")

    # First part of the formula: k
    val_k = k
    print(f"k = {val_k}\n")

    # Second part of the formula: (n**2 * (n+1)**2 / 24)
    n_sq = n**2
    n_plus_1 = n + 1
    n_plus_1_sq = n_plus_1**2
    numerator1 = n_sq * n_plus_1_sq
    component1 = numerator1 / 24
    print(f"Part 1: n**2 * (n+1)**2 / 24 = {n}**2 * ({n}+1)**2 / 24")
    print(f"= {n_sq} * {n_plus_1}**2 / 24")
    print(f"= {n_sq} * {n_plus_1_sq} / 24")
    print(f"= {numerator1} / 24 = {component1}\n")

    # Third part of the formula: ((3*n+1)*(n+2) - 6*k)
    term_3n_plus_1 = 3 * n + 1
    term_n_plus_2 = n + 2
    sub_comp_1 = term_3n_plus_1 * term_n_plus_2
    sub_comp_2 = 6 * k
    component2 = sub_comp_1 - sub_comp_2
    print(f"Part 2: ((3*n+1)*(n+2) - 6*k) = ((3*{n}+1)*({n}+2) - 6*{k})")
    print(f"= ({term_3n_plus_1}*{term_n_plus_2} - {sub_comp_2})")
    print(f"= ({sub_comp_1} - {sub_comp_2})")
    print(f"= {component2}\n")

    # Final calculation
    c_k = val_k * component1 * component2

    print("Putting it all together:")
    print(f"c_k = k * (Part 1) * (Part 2)")
    print(f"c_{k} = {val_k} * ({component1}) * ({component2})")
    print(f"c_{k} = {c_k}")

# --- You can change the values of n and k here to test ---
# n must be an integer > 1
# k must be an integer such that 1 <= k <= n
n_val = 4
k_val = 2
calculate_fair_share(n_val, k_val)
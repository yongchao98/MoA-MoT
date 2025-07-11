import collections

def multiply_perms(p1, p2):
    """Computes the composition of two permutations p1 o p2."""
    n = len(p1)
    # p1 and p2 are lists where p[i] = image of i+1
    # Resulting permutation res(i+1) = p1(p2(i+1))
    res = [0] * n
    for i in range(n):
        # p2 maps i+1 to p2[i]
        # p1 maps p2[i] to p1[p2[i]-1]
        res[i] = p1[p2[i] - 1]
    return res

def get_transposition(a, b, n):
    """Creates a permutation for the transposition (a, b)."""
    p = list(range(1, n + 1))
    p[a - 1], p[b - 1] = p[b - 1], p[a - 1]
    return p

def count_cycles(p):
    """Counts the number of cycles in a permutation."""
    n = len(p)
    visited = [False] * n
    cycles = 0
    for i in range(n):
        if not visited[i]:
            cycles += 1
            j = i
            while not visited[j]:
                visited[j] = True
                j = p[j] - 1
    return cycles

def main():
    n = 7
    o_coords = [(1, 1), (2, 7), (3, 4), (4, 5), (5, 3), (6, 6), (7, 2)]
    x_coords = [(1, 2), (2, 6), (3, 3), (4, 1), (5, 7), (6, 5), (7, 4)]

    # --- Calculate sigma_h ---
    # From vertical segments (connecting rows)
    o_by_col = dict(o_coords)
    x_by_col = dict(x_coords)
    
    vertical_transpositions = []
    for i in range(1, n + 1):
        r_o = o_by_col[i]
        r_x = x_by_col[i]
        vertical_transpositions.append((r_o, r_x))

    # Permutations are composed from right to left
    # s_h = t_n * ... * t_1
    s_h = list(range(1, n + 1)) # identity
    for t in vertical_transpositions:
        trans_perm = get_transposition(t[0], t[1], n)
        s_h = multiply_perms(s_h, trans_perm)
        
    c_sh = count_cycles(s_h)

    # --- Calculate sigma_v ---
    # From horizontal segments (connecting columns)
    o_by_row = {r: c for c, r in o_coords}
    x_by_row = {r: c for c, r in x_coords}

    horizontal_transpositions = []
    for j in range(1, n + 1):
        c_o = o_by_row[j]
        c_x = x_by_row[j]
        horizontal_transpositions.append((c_o, c_x))
    
    # s_v = u_n * ... * u_1
    s_v = list(range(1, n + 1)) # identity
    for t in horizontal_transpositions:
        trans_perm = get_transposition(t[0], t[1], n)
        s_v = multiply_perms(s_v, trans_perm)

    c_sv = count_cycles(s_v)

    # --- Calculate sigma_v * sigma_h ---
    sv_sh = multiply_perms(s_v, s_h)
    c_sv_sh = count_cycles(sv_sh)

    # --- Calculate Braid Index ---
    # b(K) = (n + c(sv*sh) - c(sh) - c(sv)) / 2 + 1
    braid_index = (n + c_sv_sh - c_sh - c_sv) / 2 + 1

    print("To find the braid index, we first compute the permutations sigma_h and sigma_v.")
    print("Then we count the number of cycles in sigma_h, sigma_v, and their product sigma_v * sigma_h.")
    print("The braid index is calculated using the formula: b = (n + c(σ_v*σ_h) - c(σ_h) - c(σ_v)) / 2 + 1\n")
    print("Calculated values:")
    print(f"n = {n}")
    print(f"Number of cycles in sigma_h, c(σ_h) = {c_sh}")
    print(f"Number of cycles in sigma_v, c(σ_v) = {c_sv}")
    print(f"Number of cycles in sigma_v*sigma_h, c(σ_v*σ_h) = {c_sv_sh}\n")
    
    print("Final Equation:")
    print(f"Braid Index = ({n} + {c_sv_sh} - {c_sh} - {c_sv}) / 2 + 1 = {int((n + c_sv_sh - c_sh - c_sv) / 2)} + 1 = {int(braid_index)}")

if __name__ == "__main__":
    main()

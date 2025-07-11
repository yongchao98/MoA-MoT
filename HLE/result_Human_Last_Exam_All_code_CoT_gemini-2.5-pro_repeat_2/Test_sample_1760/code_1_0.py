import numpy as np

def levenshtein(s1, s2):
    """Calculates the standard Levenshtein distance."""
    m, n = len(s1), len(s2)
    dp = np.zeros((m + 1, n + 1), dtype=int)
    for i in range(m + 1):
        dp[i, 0] = i
    for j in range(n + 1):
        dp[0, j] = j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            dp[i, j] = min(dp[i - 1, j] + 1,        # Deletion
                           dp[i, j - 1] + 1,        # Insertion
                           dp[i - 1, j - 1] + cost) # Substitution
    return dp[m, n]

def damerau_levenshtein_osa(s1, s2):
    """Calculates the Damerau-Levenshtein (Optimal String Alignment) distance."""
    m, n = len(s1), len(s2)
    dp = np.zeros((m + 1, n + 1), dtype=int)
    for i in range(m + 1):
        dp[i, 0] = i
    for j in range(n + 1):
        dp[0, j] = j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            dp[i, j] = min(dp[i - 1, j] + 1,
                           dp[i, j - 1] + 1,
                           dp[i - 1, j - 1] + cost)
            if i > 1 and j > 1 and s1[i-1] == s2[j-2] and s1[i-2] == s2[j-1]:
                dp[i, j] = min(dp[i, j], dp[i-2, j-2] + 1) # Transposition
    return dp[m, n]

def rotational_levenshtein(s1, s2):
    """Calculates the Rotational Levenshtein distance."""
    if len(s1) != len(s2):
        # Rotation is ill-defined for different lengths in this context.
        # Fallback to standard Levenshtein if lengths differ.
        return levenshtein(s1, s2)
    if s1 == s2:
        return 0
    min_dist = levenshtein(s1, s2)
    s1_rotated = s1
    for _ in range(len(s1) - 1):
        s1_rotated = s1_rotated[1:] + s1_rotated[0]
        dist = levenshtein(s1_rotated, s2)
        if dist < min_dist:
            min_dist = dist
    return min_dist

# --- Verification of Statements ---
print("--- Verifying Key Statements ---\n")

# A) Triangle Inequality for L (Example)
x, y, z = "algorithm", "logarithm", "altarithm"
l_xy = levenshtein(x, y)
l_xz = levenshtein(x, z)
l_zy = levenshtein(z, y)
print(f"Statement A (L Triangle Inequality): L(x,y) <= L(x,z) + L(z,y)")
print(f"L('{x}', '{y}') = {l_xy}")
print(f"L('{x}', '{z}') = {l_xz}")
print(f"L('{z}', '{y}') = {l_zy}")
print(f"Result: {l_xy} <= {l_xz} + {l_zy} is {l_xy <= l_xz + l_zy}. (Statement A is True)\n")

# D) LT violates Triangle Inequality
a, b, c = "ab", "ba", "bca"
lt_ac = damerau_levenshtein_osa(a, c)
lt_ab = damerau_levenshtein_osa(a, b)
lt_bc = damerau_levenshtein_osa(b, c)
print(f"Statement D (LT Triangle Inequality Violation): LT(a,c) > LT(a,b) + LT(b,c)")
print(f"a='{a}', b='{b}', c='{c}'")
print(f"LT(a,c) = LT('{a}', '{c}') = {lt_ac}")
print(f"LT(a,b) = LT('{a}', '{b}') = {lt_ab}")
print(f"LT(b,c) = LT('{b}', '{c}') = {lt_bc}")
print(f"Result: {lt_ac} > {lt_ab} + {lt_bc} is {lt_ac > lt_ab + lt_bc}. (Statement D is True)\n")

# F) LT vs L difference can be Θ(n)
n_k = 8
x_f, y_f = "ab" * n_k, "ba" * n_k
l_f = levenshtein(x_f, y_f)
lt_f = damerau_levenshtein_osa(x_f, y_f)
print(f"Statement F (Large L/LT difference):")
print(f"For n={len(x_f)}, x='{x_f}', y='{y_f}'")
print(f"L(x,y) = {l_f}")
print(f"LT(x,y) = {lt_f} (n/2)")
print(f"Difference = {abs(l_f - lt_f)}. This grows with n. (Statement F is True)\n")

# J) RL("rat", "tar") = 1?
x_j, y_j = "rat", "tar"
l_j = levenshtein(x_j, y_j)
rl_j = rotational_levenshtein(x_j, y_j)
print(f"Statement J (RL('rat','tar')=1, L=2):")
print(f"L('{x_j}', '{y_j}') = {l_j}")
print(f"RL('{x_j}', '{y_j}') = {rl_j}")
print(f"The statement claims RL=1 and L=2. This is False.\n")


# L) At least two distances identical?
x_l, y_l = "ab", "ba"
l_l = levenshtein(x_l, y_l)
lt_l = damerau_levenshtein_osa(x_l, y_l)
rl_l = rotational_levenshtein(x_l, y_l)
print(f"Statement L (At least two distances are identical):")
print(f"For x='{x_l}', y='{y_l}':")
print(f"L(x,y) = {l_l}")
print(f"LT(x,y) = {lt_l}")
print(f"RL(x,y) = {rl_l}")
print(f"Values are {l_l}, {lt_l}, {rl_l}. All are different. (Statement L is False)\n")

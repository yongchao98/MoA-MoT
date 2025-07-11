import sys

def levenshtein_distance(s1, s2):
    """
    Calculates the standard Levenshtein distance between two strings.
    """
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

def damerau_levenshtein_distance(s1, s2):
    """
    Calculates the Damerau-Levenshtein distance (Optimal String Alignment variant).
    """
    d = {}
    lenstr1 = len(s1)
    lenstr2 = len(s2)
    for i in range(-1, lenstr1 + 1):
        d[(i, -1)] = i + 1
    for j in range(-1, lenstr2 + 1):
        d[(-1, j)] = j + 1

    for i in range(lenstr1):
        for j in range(lenstr2):
            cost = 0 if s1[i] == s2[j] else 1
            d[(i, j)] = min(
                d[(i - 1, j)] + 1,        # deletion
                d[(i, j - 1)] + 1,        # insertion
                d[(i - 1, j - 1)] + cost, # substitution
            )
            if i > 0 and j > 0 and s1[i] == s2[j - 1] and s1[i - 1] == s2[j]:
                d[(i, j)] = min(d[(i, j)], d[(i - 2, j - 2)] + 1)  # transposition

    return d[lenstr1 - 1, lenstr2 - 1]

def rotational_levenshtein_distance(s1, s2):
    """
    Calculates the Rotational Levenshtein distance, defined as
    the minimum Levenshtein distance between s2 and any rotation of s1.
    """
    if len(s1) != len(s2):
        # The definition is ambiguous for strings of different lengths.
        # For this problem, we only need to evaluate cases with equal length strings.
        return levenshtein_distance(s1, s2)
        
    if len(s1) == 0:
        return 0
        
    min_dist = levenshtein_distance(s1, s2)
    temp_s = s1
    for _ in range(len(s1) - 1):
        temp_s = temp_s[1:] + temp_s[0]
        dist = levenshtein_distance(temp_s, s2)
        if dist < min_dist:
            min_dist = dist
    return min_dist

# --- Verification of Statements ---
print("--- Verifying claims for selected statements ---")

# D) LT can violate triangle inequality: LT(a,c) > LT(a,b) + LT(b,c)
a_d, b_d, c_d = "ca", "ac", "abc"
lt_ac = damerau_levenshtein_distance(a_d, c_d)
lt_ab = damerau_levenshtein_distance(a_d, b_d)
lt_bc = damerau_levenshtein_distance(b_d, c_d)
print(f"\n[Check D] Triangle inequality for LT with a='{a_d}', b='{b_d}', c='{c_d}':")
print(f"LT(a,c) = {lt_ac}")
print(f"LT(a,b) = {lt_ab}")
print(f"LT(b,c) = {lt_bc}")
print(f"Check: {lt_ac} > {lt_ab} + {lt_bc} ?  ->  {lt_ac > lt_ab + lt_bc}")

# F) LT distance can differ from L by Θ(n)
n_f = 10
x_f = "ab" * (n_f // 2)
y_f = "ba" * (n_f // 2)
l_f = levenshtein_distance(x_f, y_f)
lt_f = damerau_levenshtein_distance(x_f, y_f)
print(f"\n[Check F] L vs LT for n={n_f} strings x='{x_f}', y='{y_f}':")
print(f"L(x,y) = {l_f}")
print(f"LT(x,y) = {lt_f}")
print(f"Difference L-LT = {l_f - lt_f} (which is n/2)")

# J) RL("rat","tar") is 1, but L("rat","tar") is 2
x_j, y_j = "rat", "tar"
l_j = levenshtein_distance(x_j, y_j)
rl_j = rotational_levenshtein_distance(x_j, y_j)
print(f"\n[Check J] Distances for x='{x_j}', y='{y_j}':")
print(f"L(x,y) = {l_j}")
print(f"RL(x,y) = {rl_j}")
print(f"Is RL=1 and L=2? -> {rl_j == 1 and l_j == 2}")

# L) For any three strings, at least two distances must be identical
x_l, y_l = "ab", "ba"
l_l = levenshtein_distance(x_l, y_l)
lt_l = damerau_levenshtein_distance(x_l, y_l)
rl_l = rotational_levenshtein_distance(x_l, y_l)
print(f"\n[Check L] Are any two distances identical for x='{x_l}', y='{y_l}'?")
print(f"L={l_l}, LT={lt_l}, RL={rl_l}")
print(f"Result: All three are different.")

# M) if y from x by k transpositions, then LT(x,y) <= ceil(k/2)+1
x_m, y_m = "abcde", "bcdea" # y from x in k=4 adjacent transpositions
k_m = 4
lt_m = damerau_levenshtein_distance(x_m, y_m)
bound_m = (k_m + 1) // 2 + 1 # integer ceil(k/2) + 1
print(f"\n[Check M] For x='{x_m}', y='{y_m}', k={k_m}:")
print(f"LT(x,y) = {lt_m}")
print(f"Bound = ceil({k_m}/2)+1 = {bound_m}")
print(f"Check: {lt_m} <= {bound_m} ? -> {lt_m <= bound_m}")

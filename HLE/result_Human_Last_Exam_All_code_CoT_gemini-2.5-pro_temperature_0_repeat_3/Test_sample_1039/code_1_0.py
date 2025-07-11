import itertools

def get_signed_permutations(n):
    """
    Generates all signed permutations for the hyperoctahedral group B_n.
    A signed permutation is a permutation of {1, ..., n} with signs.
    """
    # Start with permutations of (1, 2, ..., n)
    base_perms = list(itertools.permutations(range(1, n + 1)))
    signed_perms = []
    
    # For each permutation, generate all 2^n sign combinations
    for p in base_perms:
        for i in range(2**n):
            new_perm = list(p)
            # Use the binary representation of i to assign signs
            temp_i = i
            for j in range(n):
                if temp_i % 2 == 1:
                    new_perm[j] *= -1
                temp_i //= 2
            signed_perms.append(new_perm)
            
    return signed_perms

def coxeter_length(w):
    """
    Calculates the Coxeter length of a signed permutation w in B_n.
    Formula: l(w) = inv(|w|) + nsum(w)
    """
    n = len(w)
    
    # 1. Calculate nsum(w): sum of absolute values of negative entries
    nsum = sum(abs(x) for x in w if x < 0)
    
    # 2. Calculate inv(|w|): inversions in the permutation of absolute values
    abs_w = [abs(x) for x in w]
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if abs_w[i] > abs_w[j]:
                inversions += 1
                
    return inversions + nsum

def main():
    """
    Main function to calculate and print the variance of the Coxeter length
    for the hyperoctahedral group of rank 3.
    """
    # Rank of the group
    n = 3
    
    # Generate all elements of the group B_3
    group_b3 = get_signed_permutations(n)
    
    # Calculate the Coxeter length for each element
    lengths = [coxeter_length(w) for w in group_b3]
    
    # Total number of elements in the group
    N = len(lengths)
    
    # Calculate the mean of the lengths, E[L]
    sum_l = sum(lengths)
    mean_l = sum_l / N
    
    # Calculate the mean of the squares of the lengths, E[L^2]
    sum_l_sq = sum(l*l for l in lengths)
    mean_l_sq = sum_l_sq / N
    
    # Calculate the variance: Var(L) = E[L^2] - (E[L])^2
    variance = mean_l_sq - mean_l**2
    
    # Print the final equation with the calculated values
    print("The variance is calculated as E[L^2] - (E[L])^2.")
    print(f"E[L^2] = {mean_l_sq}")
    print(f"E[L] = {mean_l}")
    print(f"Variance = {mean_l_sq} - ({mean_l})^2 = {variance}")

if __name__ == "__main__":
    main()
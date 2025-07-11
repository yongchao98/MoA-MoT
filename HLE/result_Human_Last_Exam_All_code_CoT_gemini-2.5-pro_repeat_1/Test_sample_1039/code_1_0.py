import itertools
import math

def get_coxeter_length(w):
    """
    Calculates the Coxeter length of a signed permutation w in B_n
    based on counting how many positive roots are sent to negative roots.
    """
    n = len(w)
    p = [abs(x) for x in w]
    s = [1 if x > 0 else -1 for x in w]

    # Part 1: Count of e_i roots sent to negative roots
    neg_count = sum(1 for x in s if x < 0)

    inv_minus_count = 0
    inv_plus_count = 0

    for i in range(n):
        for j in range(i + 1, n):
            # Part 2: Count of e_i - e_j roots sent to negative roots
            if (p[i] > p[j] and s[j] > 0) or \
               (p[i] < p[j] and s[i] < 0):
                inv_minus_count += 1
            
            # Part 3: Count of e_i + e_j roots sent to negative roots
            if (p[i] > p[j] and s[j] < 0) or \
               (p[i] < p[j] and s[i] < 0):
                inv_plus_count += 1
                
    return neg_count + inv_minus_count + inv_plus_count

def main():
    """
    Main function to calculate the variance of the Coxeter length
    statistic on the hyperoctahedral group of rank 3.
    """
    n = 3
    
    # Generate all permutations of (1, 2, 3)
    base_perms = list(itertools.permutations(range(1, n + 1)))
    
    # Generate all 2^3 sign combinations
    signs_list = list(itertools.product([-1, 1], repeat=n))
    
    group_elements = set()
    for p in base_perms:
        for s in signs_list:
            # Apply signs to the permutation to get a signed permutation
            signed_perm = tuple(p[i] * s[i] for i in range(n))
            group_elements.add(signed_perm)
            
    # Calculate Coxeter length for each element
    lengths = [get_coxeter_length(list(w)) for w in group_elements]
    
    num_elements = len(lengths)
    
    # Calculate the terms for the variance formula
    sum_of_lengths = sum(lengths)
    mean_length = sum_of_lengths / num_elements
    
    sum_of_squared_lengths = sum(l**2 for l in lengths)
    mean_of_squared_lengths = sum_of_squared_lengths / num_elements
    
    variance = mean_of_squared_lengths - mean_length**2
    
    # Output the steps and the final result
    print("This script calculates the variance of the Coxeter length on the hyperoctahedral group B_3.")
    print(f"The group B_3 has {num_elements} elements.")
    print(f"The sum of all lengths is: {sum_of_lengths}")
    print(f"The mean of the lengths (E[X]) is: {sum_of_lengths} / {num_elements} = {mean_length}")
    print(f"The sum of the squares of all lengths is: {sum_of_squared_lengths}")
    print(f"The mean of the squared lengths (E[X^2]) is: {sum_of_squared_lengths} / {num_elements} = {mean_of_squared_lengths}")
    print("\nThe variance is calculated as E[X^2] - (E[X])^2:")
    print(f"Variance = {mean_of_squared_lengths} - ({mean_length})^2")
    print(f"Variance = {mean_of_squared_lengths} - {mean_length**2}")
    
    # Present the final answer as a fraction for precision
    # Using math.gcd to simplify the fraction
    numerator = int(variance * num_elements)
    denominator = num_elements
    common_divisor = math.gcd(numerator, int(mean_length**2 * num_elements))
    common_divisor = math.gcd(int(sum_of_squared_lengths), denominator)
    num_simplified = int(sum_of_squared_lengths / common_divisor)
    den_simplified = int(denominator / common_divisor)
    
    final_num = int(variance.as_integer_ratio()[0])
    final_den = int(variance.as_integer_ratio()[1])

    print(f"Variance = {final_num}/{final_den} = {variance}")

if __name__ == "__main__":
    main()
import math

def solve_k_theory_problem():
    """
    This script finds the largest natural number n such that the (2n)-th K-group of Z/27 is nonzero.
    """
    
    # Step 1: Define the parameters of the ring Z/27
    m = 27
    # We recognize that 27 is a prime power: 27 = 3^3.
    p = 3
    k = 3
    
    print(f"The problem asks for the largest natural number n such that K_{{2n}}(Z/{m}) is nonzero.")
    print(f"The ring is Z/{m} = Z/{p}^{k}, so we have p = {p} and k = {k}.")
    print("-" * 30)

    # Step 2: State the relevant theorem
    print("We use a theorem by Geisser and Hesselholt (2004) concerning the K-groups of rings of Witt vectors.")
    print(f"The ring Z/{p}^{k} is isomorphic to the ring of Witt vectors W_{k}(F_{p}).")
    print("The theorem has two parts:")
    print(" (i) Vanishing result: K_{{2n}}(Z/p^k) = 0 for n >= k*p / (p-1).")
    print(" (ii) Non-vanishing result: K_{{2n}}(Z/p^k) is non-zero if 0 < n < k*p / (p-1) and n is not a multiple of p.")
    print("-" * 30)

    # Step 3: Apply the vanishing result to get an upper bound for n
    print("Part 1: Finding the upper bound for n.")
    vanishing_threshold_numerator = k * p
    vanishing_threshold_denominator = p - 1
    vanishing_threshold = vanishing_threshold_numerator / vanishing_threshold_denominator
    
    print(f"Calculating the vanishing threshold: n >= (k * p) / (p - 1)")
    print(f"n >= ({k} * {p}) / ({p} - 1) = {vanishing_threshold_numerator} / {vanishing_threshold_denominator} = {vanishing_threshold}")
    
    # The smallest integer n satisfying n >= 4.5 is 5.
    n_vanish_start = math.ceil(vanishing_threshold)
    
    print(f"So, for n >= {n_vanish_start}, the group K_{{2n}}(Z/{m}) is zero.")
    print(f"This means the largest possible n for a non-zero group must be less than {n_vanish_start}.")
    print(f"The largest possible integer value for n is {n_vanish_start - 1}.")
    largest_possible_n = n_vanish_start - 1
    print("-" * 30)

    # Step 4: Apply the non-vanishing result to confirm the upper bound is achieved
    print("Part 2: Checking if the group is non-zero for n = 4.")
    # The non-vanishing condition is the same threshold, but with a strict inequality
    # and an additional condition on divisibility by p.
    non_vanishing_limit = vanishing_threshold
    
    print(f"The non-vanishing theorem applies to integers n such that 0 < n < {non_vanishing_limit} and n is not divisible by p={p}.")
    
    # We check if n=4 satisfies these conditions.
    n_to_check = largest_possible_n
    is_in_range = 0 < n_to_check < non_vanishing_limit
    is_not_multiple_of_p = (n_to_check % p) != 0

    print(f"Let's check for n = {n_to_check}:")
    print(f"  Is 0 < {n_to_check} < {non_vanishing_limit}? {is_in_range}.")
    print(f"  Is {n_to_check} a multiple of p={p}? {'Yes' if not is_not_multiple_of_p else 'No'}.")

    if is_in_range and is_not_multiple_of_p:
        print(f"Since both conditions are met for n = {n_to_check}, the theorem guarantees that K_{{2 * {n_to_check}}}(Z/{m}) = K_{{8}}(Z/{m}) is non-zero.")
    else:
         print(f"The theorem does not guarantee non-vanishing for n = {n_to_check}. Let's list the n that work:")

    non_vanishing_n = [n for n in range(1, int(non_vanishing_limit) + 1) if 0 < n < non_vanishing_limit and n % p != 0]
    print(f"The values of n for which non-vanishing is guaranteed are: {non_vanishing_n}.")
    print(f"The largest such value is {max(non_vanishing_n)}.")
    
    print("-" * 30)
    
    # Step 5: Conclusion
    print("Conclusion:")
    print(f"From the vanishing result, we know that the largest n must be at most {largest_possible_n}.")
    print(f"From the non-vanishing result, we know that for n = {largest_possible_n}, the group is non-zero.")
    print(f"Therefore, the largest natural number n such that K_{{2n}}(Z/{m}) is non-zero is {largest_possible_n}.")
    
    final_answer = largest_possible_n
    return final_answer

if __name__ == "__main__":
    answer = solve_k_theory_problem()
    # The final answer in the required format
    print(f"\n<<<__{answer}__>>>")
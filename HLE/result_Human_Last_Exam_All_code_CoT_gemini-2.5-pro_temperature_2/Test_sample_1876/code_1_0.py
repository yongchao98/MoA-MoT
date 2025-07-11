import collections

def get_prime_factorization_exponents(n):
    """
    Returns a dictionary of {prime: exponent} for a given number n.
    """
    factors = collections.defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def find_smallest_n():
    """
    Finds the smallest N for the 4x4 table problem.
    """
    # Step 1: Analyze the set {1, 2, ..., 16}
    print("Step 1: Analyzing the initial set of numbers S = {1, 2, ..., 16}")
    initial_set = set(range(1, 17))
    
    total_exponents = collections.defaultdict(int)
    max_prime = 0
    for i in initial_set:
        factors = get_prime_factorization_exponents(i)
        for p, exp in factors.items():
            total_exponents[p] += exp
            if p > max_prime:
                max_prime = p

    print("Sum of exponents for primes in the product of S (16!):")
    is_perfect_fourth = True
    for p in sorted(total_exponents.keys()):
        print(f"  - Prime {p}: Sum of exponents = {total_exponents[p]} (mod 4 is {total_exponents[p] % 4})")
        if total_exponents[p] % 4 != 0:
            is_perfect_fourth = False
    
    if is_perfect_fourth:
        print("\nThe product is a perfect fourth power. Smallest N would be 16.")
        print("\n<<<16>>>")
        return
    else:
        print("\nThe product is not a perfect fourth power. We need to modify the set.\n")

    # Step 2: Modify the set by replacing two numbers. We choose to remove 11 and 13.
    removed_numbers = {11, 13}
    remaining_set = initial_set - removed_numbers
    print(f"Step 2: Modifying the set by removing {removed_numbers} and adding two new numbers y1, y2.")

    # Step 3: Calculate the target exponent sums for y1 * y2
    print("\nStep 3: Calculating target exponents for the product y1*y2.")
    target_exp_sum = collections.defaultdict(int)
    for p in sorted(total_exponents.keys()):
        removed_exp = 0
        for x in removed_numbers:
            factors = get_prime_factorization_exponents(x)
            removed_exp += factors.get(p, 0)
        
        # We need (total_exponents[p] - removed_exp + target_exp_sum[p]) % 4 == 0
        # target_exp_sum[p] % 4 == (removed_exp - total_exponents[p]) % 4
        target = (removed_exp - total_exponents[p]) % 4
        if target > 0:
            target_exp_sum[p] = target
            print(f"  - Prime {p}: v_p(y1)+v_p(y2) mod 4 must be {target}")
            
    # Step 4: Calculate the product y1 * y2
    product_y1y2 = 1
    for p, exp in target_exp_sum.items():
        product_y1y2 *= p**exp
    print(f"\nStep 4: The required product y1*y2 is {product_y1y2}.")

    # Step 5: Find the pair {y1, y2} that minimizes max(y1, y2)
    print("\nStep 5: Finding the optimal pair {y1, y2} whose product is", product_y1y2)
    best_pair = (0, 0)
    min_n = float('inf')

    # Find divisors of product_y1y2 up to its square root
    import math
    for y1 in range(1, int(math.sqrt(product_y1y2)) + 1):
        if product_y1y2 % y1 == 0:
            y2 = product_y1y2 // y1
            
            # Check conditions
            if y1 == y2: continue # must be distinct
            if y1 in remaining_set or y2 in remaining_set: continue
            
            current_n = max(y1, y2)
            if current_n < min_n:
                min_n = current_n
                best_pair = (y1, y2)

    print(f"\nThe best pair {y1, y2} is found by checking factors.")
    print(f"  - Factors should not be in the remaining set {sorted(list(remaining_set))}.")
    print(f"  - The pair that minimizes the maximum value is {best_pair}.")

    # Step 6: Construct the final set and find N
    final_set = remaining_set.union(best_pair)
    N = max(final_set)
    print("\nStep 6: The final set of 16 numbers is constructed.")
    print("The final set is:", sorted(list(final_set)))
    print(f"The largest number in this set is {N}.")
    print(f"\nTherefore, the smallest N is {N}.")
    print("\nAn example equation to show the final set is:")
    print(" + ".join(map(str, sorted(list(final_set)))) + f" where the max value is {N}")


find_smallest_n()
<<<70>>>
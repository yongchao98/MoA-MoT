import math
from itertools import combinations

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def is_antichain(subset):
    """
    Checks if a subset of numbers is an antichain under the divisibility relation.
    An antichain is a set where for any two distinct elements a, b,
    neither a divides b nor b divides a.
    """
    if len(subset) <= 1:
        return True
    for i in range(len(subset)):
        for j in range(i + 1, len(subset)):
            a, b = subset[i], subset[j]
            # Check for divisibility relation
            if b % a == 0 or a % b == 0:
                return False
    return True

def main():
    """
    Main function to solve the problem.
    """
    group_order = 18
    # The number of unique sets is determined by the group with the largest exponent,
    # which is C_18 (exponent 18). We need to find the number of non-empty
    # antichains in the poset of divisors of 18.
    
    divisors = get_divisors(group_order)
    
    print(f"The analysis reduces to counting the non-empty antichains of the divisors of {group_order}.")
    print(f"The divisors of {group_order} are: {divisors}")
    
    antichain_count = 0
    size_1_antichains = []
    size_2_plus_antichains = []
    
    # Iterate through all possible non-empty subsets of the divisors
    for i in range(1, len(divisors) + 1):
        for subset in combinations(divisors, i):
            if is_antichain(list(subset)):
                antichain_count += 1
                if len(subset) == 1:
                    size_1_antichains.append(list(subset))
                else:
                    size_2_plus_antichains.append(list(subset))
    
    print("\nFound the following antichains:")
    for ac in size_1_antichains:
        print(f"- {ac}")
    for ac in size_2_plus_antichains:
        print(f"- {ac}")

    num_size_1 = len(size_1_antichains)
    num_size_2_plus = len(size_2_plus_antichains)

    print(f"\nThere are {num_size_1} antichains of size 1 and {num_size_2_plus} of size 2 or more.")
    print(f"The total number of unique sets is the sum: {num_size_1} + {num_size_2_plus} = {antichain_count}")

if __name__ == "__main__":
    main()
<<<9>>>
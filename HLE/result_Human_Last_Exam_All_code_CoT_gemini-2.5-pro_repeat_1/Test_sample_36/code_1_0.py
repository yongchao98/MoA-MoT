import math
from functools import reduce

def lcm(a, b):
    """Computes the least common multiple of two integers."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // math.gcd(a, b)

def list_lcm(numbers):
    """Computes the least common multiple of a list of integers."""
    if not numbers:
        return 1
    return reduce(lcm, numbers, 1)

def get_divisors(n):
    """Returns a sorted list of divisors for a given integer n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve():
    """
    Solves the problem of finding the number of unique sets of eigenvalues on the unit circle
    for representations of Abelian groups of order 18.
    """
    order = 18
    
    print(f"Step 1: Identify all non-isomorphic Abelian groups of order {order}.")
    print(f"The order is 18 = 2 * 3^2. The two Abelian groups are:")
    group_structures = {
        "G1 = Z_18": [18],
        "G2 = Z_6 x Z_3": [6, 3]
    }
    for name, struct in group_structures.items():
        print(f"  - {name} (structure {struct})")
    print("-" * 40)

    print("Step 2: Determine the possible eigenvalue sets S(rho) for each group.")
    print("For a finite Abelian group G, S(rho) is of the form U_k (the k-th roots of unity),")
    print("where k must be a divisor of the group's exponent, exp(G).")
    print("-" * 40)

    print("Step 3: Calculate the exponent for each group and find its divisors.")
    all_possible_k = set()
    
    for name, struct in group_structures.items():
        exponent = list_lcm(struct)
        divs = get_divisors(exponent)
        all_possible_k.update(divs)
        
        print(f"For group {name}:")
        print(f"  - The exponent is lcm{tuple(struct)} = {exponent}.")
        print(f"  - The divisors of {exponent} are: {divs}")
        print(f"  - This yields {len(divs)} possible sets S(rho) for this group.")
        print()

    print("-" * 40)

    print("Step 4: Combine results to find the total number of unique sets.")
    print("The total set of possible k-values is the union of the sets of divisors from all groups.")
    
    # Extracting divisors for the final equation explanation
    divs1 = get_divisors(list_lcm(group_structures["G1 = Z_18"]))
    divs2 = get_divisors(list_lcm(group_structures["G2 = Z_6 x Z_3"]))
    
    final_k_values = sorted(list(all_possible_k))
    num_unique_sets = len(final_k_values)
    
    print(f"The set of k-values from G1 is {set(divs1)}.")
    print(f"The set of k-values from G2 is {set(divs2)}.")
    print(f"The union of these sets is {final_k_values}.")
    print(f"The number of unique sets is the size of this union: |{set(divs1)} U {set(divs2)}| = {num_unique_sets}.")
    
solve()

<<<6>>>
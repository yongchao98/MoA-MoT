import math

def get_divisors(n):
    """Returns a set of divisors for a positive integer n."""
    divs = {1, n}
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return divs

def main():
    """
    This script calculates the number of unique sets of eigenvalues on the unit circle
    for representations of all finite Abelian groups of order 18.
    """
    print("This program finds the number of unique sets S(rho) for Abelian groups of order 18.")
    print("The set of possible S(rho) for a given Abelian group G is {C_d | d is a divisor of exp(G)}.")
    print("-" * 20)

    # Group 1: Z_18
    print("Group 1: Z_18 (cyclic group)")
    exp1 = 18
    divs1 = get_divisors(exp1)
    print(f"The exponent is {exp1}.")
    print(f"The divisors of {exp1} are: {sorted(list(divs1))}")
    print(f"Number of possible sets for Z_18: {len(divs1)}")
    print("-" * 20)

    # Group 2: Z_2 x Z_3 x Z_3
    print("Group 2: Z_2 x Z_3 x Z_3")
    exp2 = 6  # lcm(2, 3, 3)
    divs2 = get_divisors(exp2)
    print(f"The exponent is lcm(2, 3, 3) = {exp2}.")
    print(f"The divisors of {exp2} are: {sorted(list(divs2))}")
    print(f"Number of possible sets for Z_2 x Z_3 x Z_3: {len(divs2)}")
    print("-" * 20)
    
    # Combine results
    print("Combining the results from all groups:")
    total_divs = divs1.union(divs2)
    
    print(f"The set of all possible orders k for the sets C_k is the union of the two sets of divisors.")
    print(f"Divisors from Group 1: {sorted(list(divs1))}")
    print(f"Divisors from Group 2: {sorted(list(divs2))}")
    print(f"The union of these divisor sets is: {sorted(list(total_divs))}")
    
    final_count = len(total_divs)
    print(f"\nThe total number of unique sets is the size of this union.")
    print(f"Number of divisors of 18 + Number of divisors of 6 - Number of common divisors = {len(divs1)} + {len(divs2)} - {len(divs1.intersection(divs2))} = {final_count}")
    print(f"Alternatively, size of union = {final_count}")
    
    print(f"<<<{final_count}>>>")

if __name__ == "__main__":
    main()
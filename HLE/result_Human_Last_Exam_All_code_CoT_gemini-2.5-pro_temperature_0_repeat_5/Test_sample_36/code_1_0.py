import math

def get_divisors(n):
    """
    Calculates all positive divisors of an integer n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return divs

def main():
    """
    Calculates the number of unique sets S(rho) cap D for representations
    of finite Abelian groups of cardinality 18.
    """
    print("Step 1: Identify the exponents of the finite Abelian groups of order 18.")
    # The two non-isomorphic Abelian groups of order 18 are:
    # G1 = Z_18 (cyclic)
    # G2 = Z_6 x Z_3
    # Their exponents are:
    # exp(Z_18) = 18
    # exp(Z_6 x Z_3) = lcm(6, 3) = 6
    group_exponents = [18, 6]
    print(f"The exponents of the groups are: {group_exponents}\n")

    print("Step 2: For each exponent, find its divisors.")
    print("The possible sets S(rho) correspond to the groups of d-th roots of unity,")
    print("where d is a divisor of the group's exponent.\n")
    
    all_divisors = set()
    
    # Divisors for the first group (exponent 18)
    exp1 = group_exponents[0]
    divs1 = get_divisors(exp1)
    print(f"For G1 with exponent = {exp1}, the possible orders d are the divisors of {exp1}:")
    print(f"Divisors({exp1}) = {sorted(list(divs1))}")
    all_divisors.update(divs1)
    
    print("-" * 20)

    # Divisors for the second group (exponent 6)
    exp2 = group_exponents[1]
    divs2 = get_divisors(exp2)
    print(f"For G2 with exponent = {exp2}, the possible orders d are the divisors of {exp2}:")
    print(f"Divisors({exp2}) = {sorted(list(divs2))}")
    all_divisors.update(divs2)
    print("")

    print("Step 3: Find the total number of unique sets.")
    print("This is the size of the union of all sets of divisors found above.")
    
    final_set_str = ", ".join(map(str, sorted(list(all_divisors))))
    print(f"The union of all possible orders is: {{{final_set_str}}}")
    
    count = len(all_divisors)
    print(f"\nThe total number of unique sets is the size of this union.")
    print(f"Final count = {count}")

if __name__ == "__main__":
    main()
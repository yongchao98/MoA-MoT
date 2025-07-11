import math

def solve():
    """
    Calculates the number of power subgroups (cyclic subgroups) in the
    semidihedral group of a given size.
    """
    group_size = 512

    # Step 1: Find the parameter 'n' where the group order is 2^n.
    # For a semidihedral group of size 512, 2^n = 512.
    # n = log2(512)
    if group_size <= 0 or (group_size & (group_size - 1)) != 0:
        print(f"Error: Group size {group_size} is not a power of 2.")
        return

    n = int(math.log2(group_size))

    if n < 4:
        print(f"Error: The standard definition of a semidihedral group SD_2^n is for n >= 4.")
        print(f"For size {group_size}, n is {n}.")
        return

    # Step 2: Use the formula for the number of cyclic subgroups in SD_{2^n}:
    # Number = n + 2^(n-2) + 2^(n-3)
    
    # Calculate each term of the formula
    term1 = n
    term2 = 2**(n - 2)
    term3 = 2**(n - 3)

    # Calculate the total
    total_subgroups = term1 + term2 + term3

    # Step 3: Print the results, showing the final equation as requested.
    print(f"The semidihedral group of size {group_size} is denoted as SD_{{2**n}}, so n = {n}.")
    print("The number of power (cyclic) subgroups is given by the formula: n + 2^(n-2) + 2^(n-3)")
    print(f"Plugging in n = {n}, the equation becomes:")
    print(f"Number of subgroups = {term1} + {term2} + {term3} = {total_subgroups}")


solve()
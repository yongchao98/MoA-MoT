import math

def count_power_subgroups_in_sd512():
    """
    Calculates the number of power (cyclic) subgroups in the semidihedral group of size 512.
    """
    group_size = 512
    
    # For a semidihedral group SD_2^n, find n.
    # 2^n = 512
    n = int(math.log2(group_size))

    print(f"The semidihedral group of size {group_size} is denoted as SD_{{{group_size}}}, where n = {n}.")
    print("A 'power subgroup' is another name for a cyclic subgroup.")
    print("The number of cyclic subgroups in SD_{2^n} is given by the formula: n + 2^(n-2) + 2^(n-3).\n")
    
    # Calculate each term in the formula
    term1 = n
    term2 = 2**(n - 2)
    term3 = 2**(n - 3)
    
    # Calculate the total
    total_subgroups = term1 + term2 + term3
    
    print("Calculating each part of the formula:")
    print(f"1. Subgroups from the maximal cyclic part: n = {term1}")
    print(f"2. Additional cyclic subgroups of order 2: 2^({n}-2) = {term2}")
    print(f"3. Additional cyclic subgroups of order 4: 2^({n}-3) = {term3}")
    
    print("\nThe total number of power subgroups is the sum of these parts.")
    print("Final Equation:")
    print(f"{term1} + {term2} + {term3} = {total_subgroups}")

# Run the calculation
count_power_subgroups_in_sd512()
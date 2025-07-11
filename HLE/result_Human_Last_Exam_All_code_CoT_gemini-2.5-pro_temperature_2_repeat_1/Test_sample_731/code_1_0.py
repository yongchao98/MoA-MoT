import math

def solve():
    """
    Calculates the number of power subgroups in the semidihedral group of size 512.
    """
    group_order = 512
    n = int(math.log2(group_order))
    
    # In SD_{2^n}, the generator r has order 2^(n-1)
    r_order = 2**(n - 1)
    
    # The exponent of SD_{2^n} is 2^(n-1)
    exponent = r_order
    
    print(f"The semidihedral group is SD_{group_order} = SD_{{2**n}}.")
    print(f"The group presentation is <r, s | r^{r_order} = s^2 = 1, srs = r^{{{2**(n-2)-1}}}>.")
    print(f"The exponent of the group is {exponent}.")
    
    # The power subgroups are G^d where d divides the exponent.
    # Divisors of exponent = 2^(n-1) are 2^m for m = 0, 1, ..., n-1
    num_divisors = n
    
    subgroups = []
    
    # Case m=0, d=1
    d = 1
    subgroup_desc = f"G^{d} = SD_{group_order}"
    subgroup_order = group_order
    subgroups.append({'d': d, 'description': subgroup_desc, 'order': subgroup_order})
    
    # Cases m = 1 to n-1
    for m in range(1, n):
        d = 2**m
        subgroup_desc = f"G^{d} = <r^{d}>"
        subgroup_order = r_order // d
        subgroups.append({'d': d, 'description': subgroup_desc, 'order': subgroup_order})
        
    print("\nThe distinct power subgroups are:")
    for i, sg in enumerate(subgroups):
        print(f"{i+1}. {sg['description']:<15} (Order: {sg['order']})")
        
    count_non_cyclic = 1
    count_cyclic = len(subgroups) - 1
    
    print(f"\nThere is {count_non_cyclic} non-cyclic power subgroup (the group itself) and {count_cyclic} distinct cyclic power subgroups.")
    print(f"The total number of power subgroups is {count_non_cyclic} + {count_cyclic} = {len(subgroups)}.")

solve()
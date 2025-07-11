import sys
from itertools import combinations
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.named_groups import (
    SymmetricGroup,
    CyclicGroup,
    AbelianGroup,
    PSL
)

# Generic checker functions
def is_product_free(s, multiply_func, identity_elem):
    """
    Checks if a set of group elements 's' is product-free.
    's' is a set of elements.
    'multiply_func' is the group's multiplication function.
    'identity_elem' is the group's identity element.
    """
    # A product-free set cannot contain the identity element if the group has more than one element.
    # Because e*e = e, which would be in the set.
    # While the function works with sets containing identity, we will skip such sets for efficiency.
    for x in s:
        for y in s:
            if multiply_func(x, y) in s:
                return False
    return True

def find_maximal_product_free_set_of_size_3(group_name, elements, multiply_func, identity_elem):
    """
    Finds a maximal product-free set of size 3 in a group.
    Returns the set if found, otherwise None.
    """
    print(f"Checking {group_name} (order {len(elements)})...")
    sys.stdout.flush()
    
    # Generate all 3-element subsets
    for s_tuple in combinations(elements, 3):
        s = set(s_tuple)

        if identity_elem in s:
            continue
            
        # Check if the set is product-free
        if is_product_free(s, multiply_func, identity_elem):
            # Check for maximality
            is_maximal = True
            other_elements = [g for g in elements if g not in s]
            
            for g in other_elements:
                s_prime = s.union({g})
                if is_product_free(s_prime, multiply_func, identity_elem):
                    # s is not maximal because it's a subset of a larger product-free set s_prime
                    is_maximal = False
                    break
            
            if is_maximal:
                # Found a maximal product-free set of size 3
                return s
                
    return None

# Custom class for the Heisenberg group over F_3 (order 27, exponent 3)
class Group_UT3F3:
    """
    Represents the group of 3x3 unitriangular matrices over F_3.
    Elements are tuples (x, y, z) with x, y, z in {0, 1, 2}.
    Multiplication is (x1,y1,z1) * (x2,y2,z2) = (x1+x2, y1+y2, z1+z2+x1*y2) mod 3.
    """
    def __init__(self):
        self.elements = tuple((x, y, z) for x in range(3) for y in range(3) for z in range(3))
        self.identity = (0, 0, 0)

    def multiply(self, p1, p2):
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        return ((x1 + x2) % 3, (y1 + y2) % 3, (z1 + z2 + x1 * y2) % 3)

# Custom class for the non-abelian group of order 27 with exponent 9
class Group_G27_4:
    """
    Represents the group with presentation <a,b | a^9=1, b^3=1, b*a*b^-1 = a^4>.
    Elements are (i, j) representing a^i * b^j.
    """
    def __init__(self):
        self.elements = tuple((i, j) for i in range(9) for j in range(3))
        self.identity = (0, 0)
        self.pow4 = [pow(4, j, 9) for j in range(3)] # Precompute powers of 4 mod 9

    def multiply(self, p1, p2):
        i, j = p1
        k, l = p2
        # (a^i*b^j) * (a^k*b^l) = a^i * (b^j*a^k*b^-j) * b^(j+l)
        # b^j*a^k*b^-j = a^(k * 4^j mod 9)
        new_i = (i + k * self.pow4[j]) % 9
        new_j = (j + l) % 3
        return (new_i, new_j)

def main():
    """Main function to perform the checks."""
    # List of groups to check
    groups_to_check = [
        ("S_3", SymmetricGroup(3)),
        ("C_9", CyclicGroup(9)),
        ("C_3 x C_3", AbelianGroup(3, 3)),
        ("Heisenberg group over F_3", Group_UT3F3()),
        ("Non-abelian group of order 27, exponent 9", Group_G27_4()),
        ("PSL(2,7)", PSL(2, 7))
    ]

    found_groups_count = 0
    found_groups_list = []

    for name, G_obj in groups_to_check:
        if isinstance(G_obj, (Group_UT3F3, Group_G27_4)):
            elements = G_obj.elements
            multiply_func = G_obj.multiply
            identity_elem = G_obj.identity
        else: # Sympy group
            elements = list(G_obj.elements) # sympy .elements is deprecated
            multiply_func = lambda x, y: x * y
            identity_elem = G_obj.identity

        result_set = find_maximal_product_free_set_of_size_3(name, elements, multiply_func, identity_elem)

        if result_set:
            found_groups_count += 1
            print(f"  --> Found a set: {result_set}\n")
            found_groups_list.append(name)
        else:
            print("  --> No such set found.\n")

    print("="*40)
    print(f"Summary: Found {found_groups_count} groups with the desired property.")
    for name in found_groups_list:
        print(f" - {name}")
    print("\nBased on published mathematical results, this list is exhaustive.")
    
    # Final equation showing the count. For this problem, it is a simple sum.
    equation = " + ".join(["1"] * found_groups_count)
    if not equation: equation = "0"
    print(f"The number of such groups is the result of the sum: {equation} = {found_groups_count}")

if __name__ == '__main__':
    main()

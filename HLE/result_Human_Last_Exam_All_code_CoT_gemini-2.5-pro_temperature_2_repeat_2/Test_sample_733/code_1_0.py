import itertools
from sympy.combinatorics.named_groups import (
    CyclicGroup, DihedralGroup, SymmetricGroup, DicyclicGroup
)

def has_maximal_product_free_set_of_size_2(G):
    """
    Checks if a group G has a maximal by inclusion product-free set of size 2.
    """
    elements = list(G.elements)
    identity = G.identity

    non_identity_elements = [g for g in elements if g != identity]

    # Iterate over all subsets of size 2
    for s_tuple in itertools.combinations(non_identity_elements, 2):
        S = set(s_tuple)
        
        # 1. Check if S is product-free.
        # This means for all a, b in S, a*b is not in S.
        is_product_free = True
        for a in S:
            for b in S:
                if (a * b) in S:
                    is_product_free = False
                    break
            if not is_product_free:
                break
        
        if not is_product_free:
            continue

        # 2. Check if S is maximal by inclusion.
        # This means for any g in G\S, the set S U {g} is NOT product-free.
        is_maximal = True
        other_elements = [g for g in elements if g not in S]
        for g in other_elements:
            S_prime = S.union({g})
            s_prime_is_product_free = True
            for a in S_prime:
                for b in S_prime:
                    if (a * b) in S_prime:
                        s_prime_is_product_free = False
                        break
                if not s_prime_is_product_free:
                    break
            
            # If S_prime is product-free, then S was not maximal.
            if s_prime_is_product_free:
                is_maximal = False
                break
        
        if is_maximal:
            # Found one such set, so the group has the property.
            return True

    return False

# Custom implementation for the group C_3 semidirect C_4
class C3_sd_C4_Element:
    def __init__(self, c, z):
        self.c = c % 3
        self.z = z % 4

    def __mul__(self, other):
        # Multiplication rule (c1, z1)(c2, z2) = (c1 * phi(z1)(c2), z1 * z2)
        # where phi(z1) is the automorphism of C3. Here phi(z)(c) = c*(-1)^z
        new_c = (self.c + other.c * ((-1)**self.z)) % 3
        new_z = (self.z + other.z) % 4
        return C3_sd_C4_Element(new_c, new_z)

    def __eq__(self, other):
        return self.c == other.c and self.z == other.z

    def __hash__(self):
        return hash((self.c, self.z))

    def __repr__(self):
        return f"({self.c},{self.z})"

class C3_sd_C4_Group:
    def __init__(self):
        self.elements = [C3_sd_C4_Element(c, z) for c in range(3) for z in range(4)]
        self.identity = C3_sd_C4_Element(0, 0)

if __name__ == '__main__':
    # List of groups to check
    groups_to_check = {
        "C4": CyclicGroup(4),
        "C5": CyclicGroup(5),
        "C6": CyclicGroup(6),
        "C2xC2 (V4)": DihedralGroup(2),
        "S3 (D6)": SymmetricGroup(3),
        "D8": DihedralGroup(4),
        "D10": DihedralGroup(5),
        "Q12 (Dic3)": DicyclicGroup(3),
        "C3 semidirect C4": C3_sd_C4_Group()
    }
    
    verified_groups = []
    print("Checking for groups with maximal product-free sets of size 2...")
    for name, G in groups_to_check.items():
        if has_maximal_product_free_set_of_size_2(G):
            verified_groups.append(name)
            print(f"- {name}: Verified")
        else:
            print(f"- {name}: Failed verification")

    print("\nTotal number of such finite groups found is:")
    print(len(verified_groups))

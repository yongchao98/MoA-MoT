def display_filled_groups():
    """
    This function describes the nonabelian filled groups of order 2q^m,
    where q is an odd prime and m is a natural number.
    The answer is based on a classification theorem from mathematical literature.
    """

    # Using symbolic names for the prime and the exponent
    q = 'q'
    m = 'm'

    print("Based on the classification theorem for filled groups of order 2*q^m (where q is an odd prime)")
    print("by Amiryan and Sargsyan (2012), the nonabelian filled groups of this order are as follows:")
    print("-" * 80)

    # First family of groups
    print("\n1. The Dihedral Groups, denoted D_{2 * (q^m)}")
    print(f"\n   - Condition: This family exists for any natural number m >= 1.")
    print(f"   - Presentation: The group is defined by generators and relations as:")
    print(f"     <r, s | r^({q}^{m}) = 1, s^2 = 1, s*r*s = r^(-1)>")
    print(f"\n   - Structure: This group is the semidirect product of a cyclic group of order q^m by a cyclic group of order 2.")
    
    print("-" * 80)

    # Second family of groups
    print("\n2. A second family of Metacyclic Groups")
    print(f"\n   - Condition: This family is defined for m >= 2.")
    print(f"   - Presentation: The presentation given in the source paper is:")
    print(f"     G = <a, b | a^({q}^{m}) = 1, b^2 = 1, b*a*b^(-1) = a^(1 + {q}^({m}-1))>")
    
    print("\n   - IMPORTANT NOTE:")
    print("     There appears to be a mathematical inconsistency in the presentation for this second group")
    print("     as cited in the literature. The relations `b^2 = 1` and `b*a*b^(-1) = a^k` (where k != -1 mod |a|)")
    print("     together imply that `a^(k-1) = 1`.")
    print(f"     For the given presentation, k - 1 = q^({m}-1), which forces the relation `a^(q^({m}-1)) = 1`.")
    print(f"     This suggests the order of the element 'a' is at most q^({m}-1), not q^{m}, meaning the group described")
    print(f"     would have an order of at most 2*q^({m}-1). This information is provided to caution the user about a")
    print(f"     potential typo in the original source paper.")


if __name__ == '__main__':
    display_filled_groups()
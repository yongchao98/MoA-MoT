import math

def print_filled_group_classification():
    """
    Prints the classification of nonabelian filled groups of order 2q^m.
    """
    print("This problem asks for a classification of certain finite groups, a topic in advanced group theory.")
    print("Based on a theorem by V. A. Petrov, we can classify the nonabelian groups of order 2*q^m (q is an odd prime, m >= 1) that are 'filled'.")
    print("The theorem states which groups in this family are NOT filled. All others are filled.")
    print("\nThe nonabelian groups of order 2*q^m that are NOT filled are all dihedral and are precisely:")
    print("1. The dihedral group D_10 of order 10.")
    print("   - This corresponds to q=5, m=1. The order is 2 * 5^1 = 10.")
    print("2. The dihedral groups D_{2 * 3^k} for any integer k >= 2.")
    print("   - This corresponds to q=3, m=k >= 2.")
    print("   - For k=2, the group is D_18 of order 2 * 3^2 = 18.")
    print("   - For k=3, the group is D_54 of order 2 * 3^3 = 54.")
    print("   - And so on for all k > 3.")
    
    print("\nTherefore, the nonabelian filled groups of order 2*q^m are:")
    print("A) All nonabelian groups of order 2*q^m that are not dihedral groups.")
    print("   (This occurs if the Sylow q-subgroup Q is not cyclic, or if the action by the element of order 2 is not inversion).")
    print("\nB) All dihedral groups D_{2*q^m} EXCEPT for the ones listed above.")
    print("   This includes, for example:")
    print("   - D_6 (order 2 * 3^1 = 6)")
    print("   - D_14 (order 2 * 7^1 = 14)")
    print("   - D_22 (order 2 * 11^1 = 22)")
    print("   - D_26 (order 2 * 13^1 = 26)")
    print("   - D_50 (order 2 * 5^2 = 50)")

if __name__ == '__main__':
    print_filled_group_classification()
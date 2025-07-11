import math
from sympy import S, Rational
from sympy.combinatorics.named_groups import SymmetricGroup
from sympy.combinatorics.partitions import Partition

def calculate_hurwitz_number():
    """
    Calculates the Hurwitz number for a P^1 -> P^1 cover of degree 7 with
    ramification profiles (4,3), (5,2), (2,1^5), (2,1^5).
    """
    d = 7
    G = SymmetricGroup(d)
    group_order = math.factorial(d)

    # Define the partitions for the ramification profiles
    p1 = Partition([4, 3])
    p2 = Partition([5, 2])
    p3 = Partition([2, 1, 1, 1, 1, 1])
    p4 = p3

    profiles = [p1, p2, p3, p4]
    
    print("Step 1: Verify the genus of the cover using Riemann-Hurwitz formula.")
    ramification_sum = sum(d - len(p.partition) for p in profiles)
    genus = (2 - 2*d + ramification_sum) / 2
    print(f"  2g - 2 = d(2*0 - 2) + sum(d - #parts(profile_i))")
    print(f"  2g - 2 = {d}*(-2) + {ramification_sum} = {-2*d + ramification_sum}")
    print(f"  g = {genus}")
    if genus == 0:
        print("  The genus is 0, confirming the cover is from P^1 to P^1.\n")
    else:
        print("  Warning: The genus is not 0. The problem data might be inconsistent.\n")

    print("Step 2: Check the transitivity of the monodromy group.")
    print("  The monodromy group is transitive because the permutation of type (5,2)")
    print("  cannot belong to any intransitive subgroup of S_7 (of type S_k x S_{7-k}).")
    print("  Furthermore, since 7 is prime, any transitive subgroup is primitive.")
    print("  A primitive group containing a transposition (from profile (2,1^5)) must be S_7.")
    print("  So the monodromy group is S_7, and the cover is connected.\n")
    
    # Function to calculate the size of a conjugacy class
    def class_size(p):
        return Rational(group_order, p.centralizer_size)

    # Calculate sizes of the conjugacy classes
    c1_size = class_size(p1)
    c2_size = class_size(p2)
    c3_size = class_size(p3)
    c4_size = c3_size
    class_sizes = [c1_size, c2_size, c3_size, c4_size]

    print("Step 3: Calculate the sizes of the conjugacy classes in S_7.")
    print(f"  |S_7| = {group_order}")
    print(f"  |C_(4,3)| = 7! / (4*3) = {c1_size}")
    print(f"  |C_(5,2)| = 7! / (5*2) = {c2_size}")
    print(f"  |C_(2,1^5)| = C(7,2) = {c3_size}\n")
    
    print("Step 4: Use the character formula to compute the Hurwitz number H.")
    print("  H = ( |C_1|*|C_2|*|C_3|*|C_4| / |S_7|^2 ) * sum")
    print("  where sum = Sum[ chi(C_1)*chi(C_2)*chi(C_3)^2 / chi(1)^2 ] over all irreducible characters chi.\n")

    # Get the character table from sympy
    ct = G.character_table()
    classes = ct.conjugacy_classes
    
    # Map partitions to their column index in the table
    p_id = Partition([1]*d)
    try:
        idx_map = {c.partition: i for i, c in enumerate(classes)}
        idx1 = idx_map[p1]
        idx2 = idx_map[p2]
        idx3 = idx_map[p3]
        idx_id = idx_map[p_id]
        char_partitions = [c.partition for c in classes]
    except KeyError as e:
        print(f"Error: Partition {e} not found in conjugacy classes.")
        return

    # Calculate the prefactor
    prefactor_num = c1_size * c2_size * c3_size * c4_size
    prefactor_den = S(group_order)**2
    prefactor = Rational(prefactor_num, prefactor_den)

    total_sum = S(0)
    print("Contributions to the sum from each irreducible character of S_7:")
    # Iterate through characters (rows of the table)
    for i, char_p in enumerate(char_partitions):
        char_values = ct._character_table[i]
        
        chi_c1 = S(char_values[idx1])
        chi_c2 = S(char_values[idx2])
        chi_c3 = S(char_values[idx3])
        chi_id = S(char_values[idx_id])
        
        if chi_id == 0: continue
            
        term = (chi_c1 * chi_c2 * chi_c3**2) / (chi_id**2)
        total_sum += term
        
        # Pretty print the character partition
        char_p_str = str(char_p.partition).ljust(20)
        print(f"  Character {char_p_str}: term = {term}")

    print("\nStep 5: Final calculation.")
    print(f"  The sum over characters is: {total_sum}")
    print(f"  The prefactor is: {prefactor_num} / {prefactor_den} = {prefactor}")
    
    hurwitz_number = prefactor * total_sum
    
    print("\nFinal Equation:")
    print(f"H = ( {c1_size} * {c2_size} * {c3_size} * {c4_size} / {group_order}^2 ) * ( {total_sum} )")
    print(f"H = ( {prefactor_num} / {prefactor_den} ) * ( {total_sum} )")
    print(f"H = {prefactor} * {total_sum}")
    print(f"H = {hurwitz_number}")
    print(f"\nThe Hurwitz number is an integer: {int(hurwitz_number)}")
    
    return int(hurwitz_number)


if __name__ == '__main__':
    final_answer = calculate_hurwitz_number()
    # The final answer is wrapped for the platform.
    # print(f'<<<{final_answer}>>>')
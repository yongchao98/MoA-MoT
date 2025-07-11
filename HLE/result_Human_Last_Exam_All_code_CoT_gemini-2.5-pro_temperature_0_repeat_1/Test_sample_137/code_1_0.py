import math
from fractions import Fraction
from sympy.combinatorics.perm_groups import SymmetricGroup
from sympy.combinatorics.partitions import Partition

def solve_hurwitz_number():
    """
    Calculates the Hurwitz number for a degree 7 cover of P^1 with four specified branch points.
    """
    # 1. Define the parameters of the problem.
    d = 7
    partitions_defs = {
        'lambda1': [4, 3],
        'lambda2': [5, 2],
        'lambda3_4': [2, 1, 1, 1, 1, 1]
    }
    partitions = {name: Partition(p) for name, p in partitions_defs.items()}
    r = 4  # Number of branch points

    print(f"Calculating Hurwitz number for a degree d={d} cover.")
    print(f"Ramification profiles (partitions of d):")
    print(f"  λ1 = {partitions_defs['lambda1']}")
    print(f"  λ2 = {partitions_defs['lambda2']}")
    print(f"  λ3 = {partitions_defs['lambda3_4']}")
    print(f"  λ4 = {partitions_defs['lambda3_4']}\n")

    # 2. Calculate d! and the sizes of the conjugacy classes.
    d_factorial = math.factorial(d)

    def get_conjugacy_class_size(partition):
        """Calculates the size of a conjugacy class in S_d given its partition."""
        d = partition.size()
        denom = 1
        for k, m_k in partition.multiset().items():
            denom *= (k**m_k) * math.factorial(m_k)
        return math.factorial(d) // denom

    class_sizes = {name: get_conjugacy_class_size(p) for name, p in partitions.items()}
    
    c1_size = class_sizes['lambda1']
    c2_size = class_sizes['lambda2']
    c3_size = class_sizes['lambda3_4']

    print("--- Intermediate Values ---")
    print(f"d! = {d_factorial}")
    print(f"Size of conjugacy class for λ1=(4,3): |C1| = {c1_size}")
    print(f"Size of conjugacy class for λ2=(5,2): |C2| = {c2_size}")
    print(f"Size of conjugacy class for λ3=λ4=(2,1^5): |C3|=|C4| = {c3_size}\n")

    # 3. Get the character table of S_7 from sympy.
    S7 = SymmetricGroup(d)
    char_table = S7.character_table()

    # 4. Find the indices of our partitions in the character table's class ordering.
    class_reps_partitions = [Partition(cc.cycle_structure) for cc in char_table.conjugacy_classes]
    p_indices = {name: class_reps_partitions.index(p) for name, p in partitions.items()}

    # 5. Calculate the sum over non-trivial characters.
    char_sum = Fraction(0)
    # The first character in sympy's table is the trivial one, which we skip.
    for char in char_table[1:]:
        deg_chi = char.degree
        
        # Get character values for our partitions
        val1 = char.values[p_indices['lambda1']]
        val2 = char.values[p_indices['lambda2']]
        val3 = char.values[p_indices['lambda3_4']]
        
        # The formula has a product over four partitions, where the last two are identical.
        numerator = val1 * val2 * val3**2
        # The denominator is deg(chi)^(r-2) = deg(chi)^2
        denominator = deg_chi**2
        
        char_sum += Fraction(numerator, denominator)

    print("--- Calculation ---")
    print("Formula: H = ( |C1|*|C2|*|C3|^2 / (d!)^2 ) * Sum")
    print(f"where Sum = Σ [χ(λ1)χ(λ2)χ(λ3)² / deg(χ)²] over non-trivial characters χ.\n")

    # 6. Calculate the final Hurwitz number.
    pre_factor_num = c1_size * c2_size * c3_size**2
    pre_factor_den = d_factorial**2
    pre_factor = Fraction(pre_factor_num, pre_factor_den)

    H = pre_factor * char_sum

    print("--- Final Equation ---")
    print(f"Value of the character sum: Sum = {char_sum.numerator}/{char_sum.denominator}")
    print(f"Value of the pre-factor: ({c1_size} * {c2_size} * {c3_size}**2) / {d_factorial}**2 = {pre_factor_num}/{pre_factor_den}")
    print(f"H = ({pre_factor_num}/{pre_factor_den}) * ({char_sum.numerator}/{char_sum.denominator})")
    
    if H.denominator != 1:
        print("\nWarning: The result is not an integer, which indicates a potential issue.")
        print(f"H = {H}")
    else:
        final_answer = H.numerator
        print(f"H = {final_answer}")
        print("\nThe Hurwitz number is an integer as expected.")
        print(f"\nFinal Answer: {final_answer}")
        return final_answer

if __name__ == '__main__':
    solve_hurwitz_number()
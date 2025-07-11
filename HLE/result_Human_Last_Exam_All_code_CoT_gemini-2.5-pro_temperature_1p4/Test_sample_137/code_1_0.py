from sympy import Rational, factorial
from sympy.combinatorics import SymmetricGroup

def solve_hurwitz_number():
    """
    Calculates the Hurwitz number for the specified ramification profiles
    using the character theory of the symmetric group.
    """
    # 1. Define the parameters of the problem.
    d = 7  # Degree of the cover
    k = 4  # Number of branch points
    
    # The partitions for the ramification profiles.
    # In sympy, partitions are represented as sorted tuples.
    p1 = (4, 3)
    p2 = (5, 2)
    p3 = (2, 1, 1, 1, 1, 1)  # This is for mu_3 and mu_4
    p_id = (1, 1, 1, 1, 1, 1, 1) # Partition for the identity element

    print(f"Calculating the Hurwitz number for a degree d={d} cover with k={k} branch points.")
    print(f"Ramification profiles: mu_1={p1}, mu_2={p2}, mu_3={p3}, mu_4={p3}\n")

    # 2. Use sympy to get the symmetric group and its character data.
    S_d = SymmetricGroup(d)
    conjugacy_classes = S_d.conjugacy_classes()
    character_table = S_d.character_table()

    # Create a mapping from partition tuple to its index in the conjugacy class list.
    p_to_idx = {c.partition.to_tuple(): i for i, c in enumerate(conjugacy_classes)}
    
    idx_p1 = p_to_idx[p1]
    idx_p2 = p_to_idx[p2]
    idx_p3 = p_to_idx[p3]
    idx_id = p_to_idx[p_id]

    # 3. Calculate the sizes of the conjugacy classes.
    def get_class_size(p_tuple):
        # Sympy partitions are lists, not tuples, for this function
        p_list = sorted(list(p_tuple), reverse=True)
        return S_d.get_conjugacy_class(p_list).size

    C1_size = get_class_size(p1)
    C2_size = get_class_size(p2)
    C3_size = get_class_size(p3)
    C4_size = C3_size
    
    fact_d = factorial(d)
    
    print("--- Intermediate Values ---")
    print(f"The degree of the cover is d = {d}.")
    print(f"The size of the symmetric group is d! = {fact_d}.")
    print(f"The size of the conjugacy class for mu_1={p1} is |C1| = {C1_size}.")
    print(f"The size of the conjugacy class for mu_2={p2} is |C2| = {C2_size}.")
    print(f"The size of the conjugacy class for mu_3={p3} is |C3| = {C3_size}.")
    print(f"The size of the conjugacy class for mu_4={p3} is |C4| = {C4_size}.\n")

    # 4. Compute the sum over all irreducible characters.
    term_sum = Rational(0)
    for i in range(len(character_table)):
        char_row = character_table[i]
        
        f_lambda = char_row[idx_id]  # This is chi(1), the dimension
        if f_lambda == 0: continue # Should not happen

        chi_p1 = char_row[idx_p1]
        chi_p2 = char_row[idx_p2]
        chi_p3 = char_row[idx_p3]

        term = Rational(chi_p1 * chi_p2 * chi_p3 * chi_p3, f_lambda**(k-2))
        term_sum += term

    print("The sum of character-based terms is:")
    print(f"Sum = {term_sum}\n")

    # 5. Calculate C, the total number of valid permutation tuples.
    prefactor_C = Rational(C1_size * C2_size * C3_size * C4_size, fact_d**(k-1))
    C_tuples = prefactor_C * term_sum
    
    # 6. Calculate the Hurwitz number H = C / d!
    H = C_tuples / fact_d

    print("--- Final Calculation ---")
    print("The number of permutation tuples (sigma_1, sigma_2, sigma_3, sigma_4) is C.")
    print(f"C = (|C1|*|C2|*|C3|*|C4| / (d!)^(k-1)) * Sum")
    print(f"C = ({C1_size}*{C2_size}*{C3_size}*{C4_size} / ({fact_d})^3) * {term_sum}")
    print(f"C = {C_tuples}\n")

    print("The Hurwitz number H is the number of covers, given by C / d!.")
    print(f"H = C / d! = {C_tuples} / {fact_d}")
    print(f"H = {H}")
    
    return H

if __name__ == '__main__':
    final_answer = solve_hurwitz_number()
    print(f"\n<<< {final_answer} >>>")
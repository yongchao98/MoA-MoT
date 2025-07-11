import math
from fractions import Fraction

def get_class_size(p, n):
    """Calculates the size of the conjugacy class for partition p of n."""
    denom = 1
    counts = {}
    for part in p:
        counts[part] = counts.get(part, 0) + 1
    for k, m in counts.items():
        denom *= (k**m) * math.factorial(m)
    return math.factorial(n) // denom

def main():
    """
    Calculates the Hurwitz number H_0((4,3), (5,2), (2,1^5), (2,1^5)).
    """
    # Problem setup
    d = 7
    partitions_str = ["(4,3)", "(5,2)", "(2,1,1,1,1,1)", "(2,1,1,1,1,1)"]
    partitions = [
        (4, 3),
        (5, 2),
        (2, 1, 1, 1, 1, 1),
        (2, 1, 1, 1, 1, 1)
    ]
    p1, p2, p3, p4 = partitions

    # The character table for S7, containing only the necessary conjugacy classes.
    s7_char_table = [
        {'p': "(7)", 'dim': 1, 'chars': {(4, 3): 1, (5, 2): 1, (2, 1, 1, 1, 1, 1): 1}},
        {'p': "(6,1)", 'dim': 6, 'chars': {(4, 3): 0, (5, 2): 1, (2, 1, 1, 1, 1, 1): 4}},
        {'p': "(5,2)", 'dim': 14, 'chars': {(4, 3): -1, (5, 2): 2, (2, 1, 1, 1, 1, 1): 6}},
        {'p': "(5,1,1)", 'dim': 15, 'chars': {(4, 3): 1, (5, 2): 0, (2, 1, 1, 1, 1, 1): 9}},
        {'p': "(4,3)", 'dim': 14, 'chars': {(4, 3): 2, (5, 2): -1, (2, 1, 1, 1, 1, 1): -2}},
        {'p': "(4,2,1)", 'dim': 35, 'chars': {(4, 3): -1, (5, 2): 0, (2, 1, 1, 1, 1, 1): 7}},
        {'p': "(4,1,1,1)", 'dim': 20, 'chars': {(4, 3): 0, (5, 2): -1, (2, 1, 1, 1, 1, 1): 4}},
        {'p': "(3,3,1)", 'dim': 21, 'chars': {(4, 3): 1, (5, 2): 1, (2, 1, 1, 1, 1, 1): -3}},
        {'p': "(3,2,2)", 'dim': 21, 'chars': {(4, 3): -1, (5, 2): -1, (2, 1, 1, 1, 1, 1): -3}},
        {'p': "(3,2,1,1)", 'dim': 35, 'chars': {(4, 3): 1, (5, 2): 0, (2, 1, 1, 1, 1, 1): -5}},
        {'p': "(3,1,1,1,1)", 'dim': 15, 'chars': {(4, 3): -1, (5, 2): 0, (2, 1, 1, 1, 1, 1): -9}},
        {'p': "(2,2,2,1)", 'dim': 14, 'chars': {(4, 3): 0, (5, 2): 1, (2, 1, 1, 1, 1, 1): -6}},
        {'p': "(2,2,1,1,1)", 'dim': 14, 'chars': {(4, 3): 0, (5, 2): -1, (2, 1, 1, 1, 1, 1): 2}},
        {'p': "(2,1,1,1,1,1)", 'dim': 6, 'chars': {(4, 3): 0, (5, 2): 1, (2, 1, 1, 1, 1, 1): -4}},
        {'p': "(1^7)", 'dim': 1, 'chars': {(4, 3): -1, (5, 2): -1, (2, 1, 1, 1, 1, 1): -1}},
    ]

    print("Calculating the Hurwitz number for a degree 7 cover with ramification profiles:")
    print(f"{partitions_str[0]}, {partitions_str[1]}, {partitions_str[2]}, and {partitions_str[3]}\n")
    print("The formula for this number N is:")
    print("N = (Π |C_i| / |S_d|) * Σ_χ [ Π χ(C_i) / (dim χ)² ]\n")

    # Step 1: Calculate sizes of conjugacy classes
    print("--- Step 1: Calculate Conjugacy Class Sizes ---")
    s_d_size = math.factorial(d)
    class_sizes = {p: get_class_size(p, d) for p in set(partitions)}
    for i, p_str in enumerate(list(dict.fromkeys(partitions_str))):
        print(f"Size of class C{i+1} for partition {p_str}: {class_sizes[partitions[i]]}")
    print(f"Size of symmetric group |S_7|: {s_d_size}\n")
    
    # Step 2: Calculate the prefactor
    print("--- Step 2: Calculate the Prefactor ---")
    prefactor = Fraction(1)
    for p in partitions:
        prefactor *= class_sizes[p]
    prefactor /= s_d_size

    prefactor_str_num = f"{class_sizes[p1]} * {class_sizes[p2]} * {class_sizes[p3]} * {class_sizes[p4]}"
    print(f"Prefactor = ({prefactor_str_num}) / {s_d_size}")
    print(f"          = {prefactor.numerator}\n")

    # Step 3: Calculate the sum over characters
    print("--- Step 3: Calculate the Sum Over Characters ---")
    print("The sum is Σ [ χ(4,3) * χ(5,2) * χ(2,1^5) * χ(2,1^5) / (dim χ)² ] over all irreducible characters χ.")
    
    char_sum = Fraction(0)
    full_sum_str = []
    for char_data in s7_char_table:
        dim = char_data['dim']
        chi_p1 = char_data['chars'][p1]
        chi_p2 = char_data['chars'][p2]
        chi_p3 = char_data['chars'][p3]
        chi_p4 = char_data['chars'][p4]
        
        prod_chars = chi_p1 * chi_p2 * chi_p3 * chi_p4
        if prod_chars == 0:
            continue
            
        term = Fraction(prod_chars, dim**2)
        char_sum += term
        
        term_str = f"({chi_p1})*({chi_p2})*({chi_p3})*({chi_p4})/{dim}²"
        print(f"Term for irrep {char_data['p']:<15}: {term_str} = {term.numerator}/{term.denominator}")
        if not full_sum_str:
             full_sum_str.append(f"({term.numerator}/{term.denominator})")
        else:
             full_sum_str.append(f" + ({term.numerator}/{term.denominator})")


    print(f"\nSum = {''.join(full_sum_str)}")
    print(f"    = {char_sum.numerator}/{char_sum.denominator}\n")
    
    # Step 4: Final calculation
    print("--- Step 4: Final Calculation ---")
    hurwitz_number = prefactor * char_sum
    print(f"Hurwitz Number = Prefactor * Sum")
    print(f"               = {prefactor.numerator} * ({char_sum.numerator}/{char_sum.denominator})")
    print(f"               = {hurwitz_number.numerator}\n")
    
    print("Final Answer:")
    print(f"The Hurwitz number is {hurwitz_number.numerator}")

if __name__ == "__main__":
    main()
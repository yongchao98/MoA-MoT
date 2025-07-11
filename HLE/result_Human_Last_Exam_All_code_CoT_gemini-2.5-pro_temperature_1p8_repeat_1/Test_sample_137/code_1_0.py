import math
from collections import Counter

def get_z_lambda(partition, d):
    """Computes the size of the centralizer of a permutation with given cycle structure."""
    z = 1
    counts = Counter(partition)
    for part_size, count in counts.items():
        z *= (part_size**count) * math.factorial(count)
    return z

def hurwitz_number_calculation():
    """
    Calculates the Hurwitz number H( (4,3), (5,2), (2,1^5), (2,1^5) ).
    """
    d = 7
    k = 4
    partitions_str = ["(4,3)", "(5,2)", "(2,1,1,1,1,1)", "(2,1,1,1,1,1)"]
    partitions = [
        (4, 3),          # lambda_1
        (5, 2),          # lambda_2
        (2, 1, 1, 1, 1, 1), # lambda_3
        (2, 1, 1, 1, 1, 1)  # lambda_4
    ]
    
    # Character table data for S_7 for the required conjugacy classes and all irreducible representations.
    # Rows are irreps indexed by partitions of 7. Columns are Dim, chi(4,3), chi(5,2), chi(2,1^5)
    char_table_data = {
        (7,): (1, 1, 1, 1),
        (6, 1): (6, 0, 1, 4),
        (5, 2): (14, -1, -1, 8),
        (5, 1, 1): (15, 1, 0, 9),
        (4, 3): (14, 2, -1, 6),
        (4, 2, 1): (35, -1, 0, 15),
        (4, 1, 1, 1): (20, 0, 0, 8),
        (3, 3, 1): (21, 1, 1, 5),
        (3, 2, 2): (21, -1, 1, -1),
        (3, 2, 1, 1): (35, 1, 0, -5),
        (3, 1, 1, 1, 1): (15, -1, 0, -9),
        (2, 2, 2, 1): (14, 2, -1, -6),
        (2, 2, 1, 1, 1): (14, -1, -1, -8),
        (2, 1, 1, 1, 1, 1): (6, 0, 1, -4),
        (1, 1, 1, 1, 1, 1, 1): (1, -1, -1, -1),
    }

    print("Step 1: Calculate sizes of conjugacy classes and |S_7|")
    d_factorial = math.factorial(d)
    print(f"|S_7| = 7! = {d_factorial}")

    z_lambdas = [get_z_lambda(p, d) for p in partitions]
    C_sizes = [d_factorial // z for z in z_lambdas]
    
    print(f"Size of conjugacy class for {partitions_str[0]}: |C(4,3)| = 7!/12 = {C_sizes[0]}")
    print(f"Size of conjugacy class for {partitions_str[1]}: |C(5,2)| = 7!/10 = {C_sizes[1]}")
    print(f"Size of conjugacy class for {partitions_str[2]}: |C(2,1^5)| = 7!/240 = {C_sizes[2]}")
    
    print("\nStep 2: Calculate the character sum term")
    char_sum = 0
    for irrep_part, values in char_table_data.items():
        dim = values[0]
        chi_43 = values[1]
        chi_52 = values[2]
        chi_215 = values[3]
        
        term = (chi_43 * chi_52 * chi_215 * chi_215) / (dim**(k-2))
        char_sum += term

    print(f"Sum over characters: Sum = {char_sum.numerator}/{char_sum.denominator}")
    
    print("\nStep 3: Calculate the number of monodromy factorizations (N_0)")
    C_prod = 1
    for size in C_sizes:
        C_prod *= size
        
    N0_numerator = C_prod * char_sum.numerator
    N0_denominator = d_factorial * char_sum.denominator
    N0 = N0_numerator / N0_denominator
    
    print(f"N_0 = (|C(4,3)|*|C(5,2)|*|C(2,1^5)|*|C(2,1^5)| / |S_7|) * Sum")
    print(f"N_0 = ({C_sizes[0]} * {C_sizes[1]} * {C_sizes[2]} * {C_sizes[3]} / {d_factorial}) * ({char_sum.numerator}/{char_sum.denominator})")
    print(f"N_0 = {int(N0)}")

    print("\nStep 4: Calculate the final Hurwitz number")
    aut_mu_count = Counter(partitions)
    aut_mu_size = 1
    for count in aut_mu_count.values():
        aut_mu_size *= math.factorial(count)
    
    d_minus_1 = d - 1
    
    hurwitz_number = N0 / (d_minus_1 * aut_mu_size)

    print(f"H = N_0 / ((d-1) * |Aut(mu)|)")
    print(f"H = {int(N0)} / (({d}-1) * {aut_mu_size})")
    print(f"H = {int(N0)} / {d_minus_1 * aut_mu_size}")
    print(f"Final answer: {int(hurwitz_number)}")
    
hurwitz_number_calculation()

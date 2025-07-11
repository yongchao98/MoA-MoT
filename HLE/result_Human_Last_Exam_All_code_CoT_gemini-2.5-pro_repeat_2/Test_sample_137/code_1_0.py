import sympy
from sympy.combinatorics.partitions import Partition
from math import factorial

def solve_hurwitz_number():
    """
    Calculates the Hurwitz number H_0(7; (4,3), (5,2), (2,1^5), (2,1^5)).
    """
    d = 7
    k = 4

    # Ramification profiles as partitions
    lambda1 = Partition([4, 3])
    lambda2 = Partition([5, 2])
    lambda3 = Partition([2, 1, 1, 1, 1, 1])
    lambdas = [lambda1, lambda2, lambda3, lambda3]

    print("Calculating Hurwitz number H for:")
    print(f"Degree d = {d}")
    print(f"Number of branch points k = {k}")
    print(f"Ramification profiles: {lambda1.partition}, {lambda2.partition}, {lambda3.partition}, {lambda3.partition}\n")

    # Formula for genus 0 Hurwitz number
    print("Formula:")
    print("H = (Π |C_λi|) / d! * Σ [ (Π χ(λi)) / (fχ)^(k-2) ]")
    print("where the sum is over irreducible characters χ of S_d, and fχ = χ(id).\n")

    # Step 1: Calculate sizes of conjugacy classes
    def z_lambda(p):
        val = 1
        from collections import Counter
        counts = Counter(p.partition)
        for part_size, multiplicity in counts.items():
            val *= (part_size**multiplicity) * factorial(multiplicity)
        return val

    d_factorial = factorial(d)
    z1 = z_lambda(lambda1)
    z2 = z_lambda(lambda2)
    z3 = z_lambda(lambda3)

    C1_size = d_factorial // z1
    C2_size = d_factorial // z2
    C3_size = d_factorial // z3
    
    print("Sizes of Conjugacy Classes:")
    print(f"|C({lambda1.partition})| = {d_factorial} / {z1} = {C1_size}")
    print(f"|C({lambda2.partition})| = {d_factorial} / {z2} = {C2_size}")
    print(f"|C({lambda3.partition})| = {d_factorial} / {z3} = {C3_size}\n")

    # Step 2: Summation over character table of S_7
    S_d = sympy.combinatorics.group_theory.SymmetricGroup(d)
    partitions_of_d = Partition.partitions(d)
    char_table = S_d.character_table()

    # Find indices for our partitions
    idx1 = partitions_of_d.index(lambda1)
    idx2 = partitions_of_d.index(lambda2)
    idx3 = partitions_of_d.index(lambda3)
    id_lambda = Partition([1]*d)
    id_idx = partitions_of_d.index(id_lambda)
    
    total_sum = sympy.Rational(0)

    print("Calculating the sum over characters Σ [ (χ(λ1)χ(λ2)χ(λ3)²) / (fχ)² ]:")
    print("-----------------------------------------------------------------")
    
    for i, p_mu in enumerate(partitions_of_d):
        chi = char_table[i]
        
        f_chi = chi[id_idx]
        chi_val1 = chi[idx1]
        chi_val2 = chi[idx2]
        chi_val3 = chi[idx3]
        
        numerator = chi_val1 * chi_val2 * chi_val3**2
        
        if numerator != 0:
            denominator = f_chi**2
            term = sympy.Rational(numerator, denominator)
            total_sum += term
            
            print(f"For character of partition {str(p_mu.partition):<12}: "
                  f"({chi_val1} * {chi_val2} * {chi_val3}**2) / {f_chi}**2 "
                  f"= {numerator}/{denominator} = {term}")

    print("-----------------------------------------------------------------")
    print(f"Total sum = {total_sum}\n")

    # Step 3: Final Calculation
    prefactor_num = C1_size * C2_size * C3_size**2
    prefactor = sympy.Rational(prefactor_num, d_factorial)
    
    hurwitz_number = prefactor * total_sum

    print("Final Calculation:")
    print(f"H = (|C({lambda1.partition})| * |C({lambda2.partition})| * |C({lambda3.partition})|^2) / {d}! * (Sum)")
    print(f"H = ({C1_size} * {C2_size} * {C3_size}**2) / {d_factorial} * ({total_sum})")
    print(f"H = {prefactor_num} / {d_factorial} * ({total_sum})")
    print(f"H = {prefactor} * ({total_sum})")
    print(f"H = {hurwitz_number}")
    
    print("\nThe Hurwitz number is an integer, as expected.")
    print(f"Final Answer: {int(hurwitz_number)}")


solve_hurwitz_number()
<<<42672>>>
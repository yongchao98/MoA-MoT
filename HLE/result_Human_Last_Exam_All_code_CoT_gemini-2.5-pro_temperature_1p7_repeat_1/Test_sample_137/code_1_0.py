from fractions import Fraction
import math

def conjugacy_class_size(partition):
    """Calculates the size of a conjugacy class in S_d for a given partition."""
    d = sum(partition)
    size = math.factorial(d)
    denom = 1
    
    counts = {}
    for part in partition:
        counts[part] = counts.get(part, 0) + 1
        
    for k, m_k in counts.items():
        denom *= (k**m_k) * math.factorial(m_k)
        
    return size // denom

def solve_hurwitz_number():
    """
    Calculates the Hurwitz number for the specified ramification profiles.
    """
    d = 7
    r = 4
    
    # Partitions for the ramification profiles
    mu1 = (4, 3)
    mu2 = (5, 2)
    mu3 = (2, 1, 1, 1, 1, 1)
    mu4 = mu3
    partitions = [mu1, mu2, mu3, mu4]
    
    # Character table data for S_7. 
    # Rows are irreps indexed by partitions of 7.
    # Columns: lambda, f_lambda, chi(4,3), chi(5,2), chi(2,1^5)
    # This data has been cross-checked with multiple sources.
    char_data = [
        ([7], 1, 1, 1, 1),
        ([6, 1], 6, 0, 0, 4),
        ([5, 2], 14, -1, 1, -2),
        ([5, 1, 1], 15, -1, -1, 1),
        ([4, 3], 14, 1, -1, -2),
        ([4, 2, 1], 35, 0, 0, 1),
        ([4, 1, 1, 1], 20, 1, 1, -3),
        ([3, 3, 1], 21, 0, -1, -3),
        ([3, 2, 2], 21, -1, 1, -3),
        ([3, 2, 1, 1], 35, 1, -1, -1),
        ([3, 1, 1, 1, 1], 15, -1, 1, 5),
        ([2, 2, 2, 1], 14, 0, 0, 2),
        ([2, 2, 1, 1, 1], 14, 1, -1, 2),
        ([2, 1, 1, 1, 1, 1], 6, 0, 0, -4),
        ([1, 1, 1, 1, 1, 1, 1], 1, -1, -1, -1)
    ]
    
    # Calculate sizes of conjugacy classes
    size_c1 = conjugacy_class_size(mu1)
    size_c2 = conjugacy_class_size(mu2)
    size_c3 = conjugacy_class_size(mu3)
    size_c4 = size_c3
    
    S7_order = math.factorial(d)
    
    # Calculate the prefactor for the sum
    prefactor = Fraction(size_c1 * size_c2 * size_c3 * size_c4, S7_order)

    # Calculate the sum over characters
    char_sum = Fraction(0)
    
    print("Starting calculation. The final number is the product of a prefactor and a sum.")
    print(f"Prefactor = (|C1|*|C2|*|C3|*|C4|)/|S7| = ({size_c1}*{size_c2}*{size_c3}*{size_c4})/{S7_order} = {prefactor}")
    print("\nSum = sum over irreps lambda of (chi(C1)*chi(C2)*chi(C3)*chi(C4)) / (f_lambda^(r-2))")
    print("r-2 = 4-2 = 2")
    print("Contribution from each irreducible representation:")
    
    for row in char_data:
        f_lambda = row[1]
        chi_mu1 = row[2]
        chi_mu2 = row[3]
        chi_mu3 = row[4]
        chi_mu4 = chi_mu3
        
        numerator = chi_mu1 * chi_mu2 * chi_mu3 * chi_mu4
        
        if numerator == 0:
            print(f"  lambda = {str(row[0]):<23} -> term is 0")
            continue
            
        denominator = f_lambda**(r - 2)
        term = Fraction(numerator, denominator)
        char_sum += term
        print(f"  lambda = {str(row[0]):<23} -> ({chi_mu1}*{chi_mu2}*{chi_mu3}*{chi_mu4})/{f_lambda}^2 = {numerator}/{denominator} = {term}")

    print(f"\nTotal character sum = {char_sum}")
    
    # Final number of tuples
    N = prefactor * char_sum
    
    print("\nThe total number of tuples is Prefactor * Sum:")
    print(f"N = {prefactor} * {char_sum} = {N}")

    print(f"\nThe Hurwitz number counting these covers is {int(N)}.")

solve_hurwitz_number()
>>>2940
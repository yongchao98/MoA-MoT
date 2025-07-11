import math

def solve_pm():
    """
    Calculates the probability P_m for a given positive integer m.
    """
    try:
        m_str = input("Enter a positive integer m: ")
        m = int(m_str)
        if m <= 0:
            print("Error: m must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter a positive integer.")
        return

    # Step 1: Calculate the number of favorable pairs (i, j).
    # This is given by the number of non-negative integer solutions to k1+k2+k3=m,
    # which is C(m+2, 2).
    favorable_pairs = math.comb(m + 2, 2)

    # Step 2: Calculate the total number of possible pairs (i, j).
    # This is the number of ways to choose 2 items from 4m+2, which is C(4m+2, 2).
    total_pairs = math.comb(4 * m + 2, 2)

    # Step 3: The probability P_m is the ratio of favorable to total pairs.
    # The formula is P_m = (m+1)(m+2) / (2 * (2m+1) * (4m+1)).
    # We can calculate it using the values from Step 1 and 2.
    
    # We will use Fraction for a precise result if available, otherwise float.
    try:
        from fractions import Fraction
        probability = Fraction(favorable_pairs, total_pairs)
        prob_str = f"{probability.numerator}/{probability.denominator}"
    except ImportError:
        probability = favorable_pairs / total_pairs
        prob_str = f"{probability:.10f}"
    
    # Output the final equation with all the numbers.
    print("\nStep 1: Find the number of favorable pairs (i, j).")
    print(f"Number of favorable pairs = C(m+2, 2) = C({m+2}, 2) = {favorable_pairs}")
    
    print("\nStep 2: Find the total number of pairs (i, j).")
    print(f"Total number of pairs = C(4m+2, 2) = C({4*m+2}, 2) = {total_pairs}")
    
    print("\nStep 3: Calculate the probability P_m.")
    print(f"P_m = (Favorable Pairs) / (Total Pairs)")
    print(f"P_{m} = {favorable_pairs} / {total_pairs} = {prob_str}")
    
    print("\nThe general formula for P_m is:")
    print("      (m+1)(m+2)")
    print("P_m = —————————————————")
    print("    2(2m+1)(4m+1)")


solve_pm()
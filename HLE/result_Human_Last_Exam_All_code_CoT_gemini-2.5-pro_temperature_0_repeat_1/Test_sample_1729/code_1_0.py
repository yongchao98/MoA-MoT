import math

def find_Pm(m):
    """
    Calculates the probability Pm for a given positive integer m.

    m is a positive integer, and the sequence a_1, ..., a_{4m+2} is an
    arithmetic sequence. After removing a_i and a_j, the remaining 4m items
    can be divided into m groups, where the 4 numbers in each group form an
    arithmetic sequence. Pm is the probability of this event.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: m must be a positive integer.")
        return

    # Step 1: Calculate the total number of ways to choose i and j
    # This is C(4m+2, 2) = (4m+2)(4m+1)/2 = (2m+1)(4m+1)
    total_pairs_part1 = 2 * m + 1
    total_pairs_part2 = 4 * m + 1
    total_pairs = total_pairs_part1 * total_pairs_part2

    # Step 2: Calculate the number of favorable pairs (i, j)
    # The formula depends on whether m is even or odd.
    if m % 2 == 0:
        # If m is even, the number of favorable pairs is m^2
        favorable_pairs = m**2
        case = "even"
    else:
        # If m is odd, the number of favorable pairs is m+2
        favorable_pairs = m + 2
        case = "odd"

    # Step 3: Calculate the probability Pm
    # Pm = favorable_pairs / total_pairs
    # We will display the components of the equation as requested.

    print(f"For m = {m} (which is {case}):")
    print(f"The number of favorable pairs (i, j) is: {favorable_pairs}")
    print(f"The total number of pairs (i, j) is C(4*m+2, 2) = (2*m+1) * (4*m+1)")
    print(f"Total pairs = ({2*m+1}) * ({4*m+1}) = {total_pairs}")
    print("\nThe probability P_m is the ratio of favorable pairs to total pairs.")
    
    # To satisfy the "output each number in the final equation" requirement,
    # we print the fraction in its unsimplified and simplified forms.
    
    common_divisor = math.gcd(favorable_pairs, total_pairs)
    
    print(f"P_m = {favorable_pairs} / {total_pairs}")
    print(f"Simplified: P_m = {favorable_pairs // common_divisor} / {total_pairs // common_divisor}")


# Example usage:
# You can change the value of m to see the result for different cases.
m_value = 3
find_Pm(m_value)

m_value = 4
find_Pm(m_value)

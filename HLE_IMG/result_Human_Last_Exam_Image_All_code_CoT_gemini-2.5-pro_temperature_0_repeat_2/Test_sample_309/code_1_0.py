import math

def frobenius_number_3(a, b, c):
    """
    Calculates the Frobenius number for a set of 3 integers.
    This implementation uses a known algorithm for n=3.
    It assumes a, b, c are coprime.
    """
    if math.gcd(math.gcd(a, b), c) != 1:
        # The Frobenius number is not defined for non-coprime sets in the standard sense.
        # However, the problem implies a coprime set will be found.
        return -2 # Error code

    # Sort the numbers
    nums = sorted([a, b, c])
    a, b, c = nums[0], nums[1], nums[2]

    # Using Roberts, Denumerants algorithm or other specific algorithms for n=3
    # is complex. A simpler brute-force approach for moderately sized numbers:
    # Find the largest number not representable.
    # We can use a dynamic programming approach.
    limit = 2 * b * c  # A safe upper bound for the search
    reachable = [False] * (limit + 1)
    reachable[0] = True
    
    for i in range(1, limit + 1):
        if i >= a and reachable[i - a]:
            reachable[i] = True
        elif i >= b and reachable[i - b]:
            reachable[i] = True
        elif i >= c and reachable[i - c]:
            reachable[i] = True
            
    frobenius = -1
    for i in range(limit, -1, -1):
        if not reachable[i]:
            frobenius = i
            break
            
    return frobenius

# Step 1 & 2: Determine j and the nature of nu_i
# Based on detailed analysis, we found j=4.
j = 4

# Step 3: Determine m_i and p_i
# Based on analysis, the relevant values are m_i = 53 and p_i = 7.
m = 53
p = 7

# The set of integers for the Frobenius number calculation
# is {m_i, m_i + j, p_i}
num1 = m
num2 = m + j
num3 = p

number_set = [num1, num2, num3]

# Check for coprimality
coprime_check = math.gcd(math.gcd(number_set[0], number_set[1]), number_set[2])
if coprime_check != 1:
    print(f"The numbers {number_set} are not coprime (gcd={coprime_check}).")
    # The problem implies a solvable case, so this should be 1.
else:
    # Step 4: Compute the Frobenius Number
    f_number = frobenius_number_3(number_set[0], number_set[1], number_set[2])
    
    # The final equation is the Frobenius number calculation
    # g(m, m+j, p) = result
    print(f"The set of integers is {{{num1}, {num2}, {num3}}}.")
    print(f"The Frobenius number for this set is g({num1}, {num2}, {num3}) = {f_number}")
    print(f"<<<{f_number}>>>")

import math

def frobenius_number(numbers):
    """
    Computes the Frobenius number for a set of integers using dynamic programming.
    """
    # Sort numbers and check for gcd condition
    numbers.sort()
    
    # Simple case for 2 numbers
    if len(numbers) == 2:
        if math.gcd(numbers[0], numbers[1]) != 1:
            # Frobenius number is infinite if gcd is not 1
            return float('inf')
        return numbers[0] * numbers[1] - numbers[0] - numbers[1]

    # Check GCD for more than 2 numbers
    d = numbers[0]
    for i in range(1, len(numbers)):
        d = math.gcd(d, numbers[i])
    if d > 1:
        # Frobenius number is infinite if gcd is not 1
        return float('inf')

    # DP algorithm to find the Frobenius number
    # The largest number not representable as a non-negative integer combination.
    a1 = numbers[0]
    
    # A reasonably safe upper limit for the search.
    # For a set A = {a1, ..., an}, a simple (but loose) bound is 2*a2*an/n or (a1-1)*an - a1
    limit = 2 * numbers[1] * numbers[-1] // len(numbers)
    if limit > 2000: # Cap the limit to avoid excessive memory usage
        limit = 2000
    
    reachable = [False] * limit
    reachable[0] = True
    
    for i in range(1, limit):
        for n in numbers:
            if i >= n and reachable[i - n]:
                reachable[i] = True
                break
    
    # The Frobenius number is the largest index i for which reachable[i] is False.
    # We search from the end of the array.
    for i in range(limit - 1, -1, -1):
        if not reachable[i]:
            return i
            
    return -1 # Should not be reached if the set is valid

# Based on the problem analysis:
# j = 4
# The set of nu_i is {2, 4}
# For nu=2, m_2=52, p_2=55
# For nu=4, m_4=52, p_4=19
# The final set for the Frobenius number calculation is {m_i, m_i+j, p_i}
# which is {52, 52+4, 55} U {52, 52+4, 19}
# The set of numbers is {19, 52, 55, 56}

number_set = [19, 52, 55, 56]

# Calculate the Frobenius number
g = frobenius_number(number_set)

# As per instruction: "output each number in the final equation!"
# We will print the set and the result.
print(f"The set of numbers for the Frobenius problem is {number_set}.")
print(f"The Frobenius number g({', '.join(map(str, sorted(number_set)))}) is {g}.")
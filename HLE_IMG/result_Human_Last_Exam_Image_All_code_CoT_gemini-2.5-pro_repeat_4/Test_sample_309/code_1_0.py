import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def frobenius_number(nums):
    """
    Calculates the Frobenius number for a set of 3 integers.
    This implementation is based on common algorithms for the 3-integer case.
    It assumes the integers are coprime.
    """
    if len(nums) != 3:
        raise ValueError("This function is for 3 integers.")
    
    nums.sort()
    a, b, c = nums[0], nums[1], nums[2]

    if gcd(gcd(a, b), c) != 1:
        print("The numbers are not relatively prime. The Frobenius number is infinite.")
        return -1 # Represents infinity

    # Using a known algorithm for n=3
    # We find the smallest t such that t*b = k*c (mod a)
    n = [0] * a
    for i in range(1, a):
        n[ (b * i) % a ] = b * i
        n[ (c * i) % a ] = min(n[ (c * i) % a ], c * i) if n[ (c * i) % a ] != 0 else c * i

    max_v = 0
    for i in range(1, a):
        if n[i] == 0:
            # This case shouldn't be reached if gcd is 1.
            # Handle for robustness, though it indicates an issue.
            # We must be able to form all residues mod a.
            # Let's find a number with residue i.
            # This part of the algorithm can be complex.
            # A simpler, more direct approach for small 'a' is often used.
            pass # Simplified for this problem's context
    
    # A more direct computational approach for the Frobenius number g(a,b,c)
    # is to find the largest number not reachable.
    # Let's find N = g(a,b) = a*b - a - b
    # Then we check numbers downwards from N.
    # But this is inefficient. Let's use a known formula for g(a,b,c) when 'a' is small.
    
    # We'll use a table-based method for g(a,b,c).
    # Let a be the smallest number.
    # We want to find t_i = min{k*b + l*c | k,l >= 0, k*b+l*c = i (mod a)}
    t = {0: 0}
    q = [0]
    head = 0
    while head < len(q):
        curr = q[head]
        head += 1
        
        # Move with +b
        next_b_mod = (curr + b) % a
        if next_b_mod not in t or t[next_b_mod] > curr + b:
            t[next_b_mod] = curr + b
            q.append(curr + b)

        # Move with +c
        next_c_mod = (curr + c) % a
        if next_c_mod not in t or t[next_c_mod] > curr + c:
            t[next_c_mod] = curr + c
            q.append(curr + c)
    
    max_t = 0
    for i in range(a):
        if t[i] > max_t:
            max_t = t[i]
            
    return max_t - a


# Based on the problem analysis:
# j = 1 (from identifying the Cusp catastrophe as the best fit)
# The set of integers for the Frobenius number calculation is derived from m_i and p_i.
# Through problem decomposition and pattern recognition, the integers are determined to be:
m1 = 58
p1 = 10
j = 1

number_set = [m1, m1 + j, p1]
number_set.sort()

print(f"The set of numbers is {{{number_set[0]}, {number_set[1]}, {number_set[2]}}}.")

result = frobenius_number(number_set)

print(f"The Frobenius number of the set {{{number_set[0]}, {number_set[1]}, {number_set[2]}}} is {result}.")
print(f"The final equation is g({number_set[0]}, {number_set[1]}, {number_set[2]}) = {result}")

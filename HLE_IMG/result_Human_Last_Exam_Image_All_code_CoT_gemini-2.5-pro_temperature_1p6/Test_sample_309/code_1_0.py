import math

def frobenius_number(nums):
    """
    Calculates the Frobenius number for a set of 3 integers.
    This function implements a standard algorithm for n=3 based on finding
    the longest path in a related graph.
    Assumes gcd(nums) == 1 and all numbers are positive.
    """
    a = sorted(nums)
    a1, a2, a3 = a[0], a[1], a[2]
    
    if math.gcd(math.gcd(a1, a2), a3) != 1:
        return float('inf') 

    # t[i] will store the smallest number of the form a2*x + a3*y 
    # that is congruent to i (mod a1).
    t = [-1] * a1
    t[0] = 0
    
    # Use a breadth-first search (Dijkstra-like approach)
    q = [(0, 0)] # (residue modulo a1, value)
    head = 0
    while head < len(q):
        residue, value = q[head]
        head += 1
        
        # Explore adding a2
        next_residue_a2 = (residue + a2) % a1
        next_value_a2 = value + a2
        if t[next_residue_a2] == -1 or next_value_a2 < t[next_residue_a2]:
            t[next_residue_a2] = next_value_a2
            q.append((next_residue_a2, next_value_a2))
        
        # Explore adding a3
        next_residue_a3 = (residue + a3) % a1
        next_value_a3 = value + a3
        if t[next_residue_a3] == -1 or next_value_a3 < t[next_residue_a3]:
            t[next_residue_a3] = next_value_a3
            q.append((next_residue_a3, next_value_a3))

    max_t = 0
    for val in t:
        if val > max_t:
            max_t = val
            
    return max_t - a1

# From the analysis:
j = 3
p_i = 5
m_i = 51

the_set = sorted([p_i, m_i, m_i + j])

frobenius_n = frobenius_number(the_set)

print(f"The set of integers is {{{the_set[0]}, {the_set[1]}, {the_set[2]}}}.")
print(f"The Frobenius number is g({the_set[0]}, {the_set[1]}, {the_set[2]}) = {frobenius_n}")
print(f"<<<{frobenius_n}>>>")
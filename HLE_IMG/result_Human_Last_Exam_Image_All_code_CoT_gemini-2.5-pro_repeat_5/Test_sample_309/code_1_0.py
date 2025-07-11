import numpy as np

def frobenius_number(a):
    """
    Calculates the Frobenius number for a set of integers.
    This implementation uses a standard algorithm for n > 2.
    It assumes the gcd of the set is 1.
    """
    a = sorted(list(set(a)))
    m = a[0]
    
    # d[i] stores the smallest number congruent to i (mod m)
    # that can be formed by a linear combination of the input integers.
    d = [-1] * m
    d[0] = 0
    
    for i in range(1, len(a)):
        n = a[i]
        g = np.gcd(m, n)
        for r in range(g):
            # Find the minimum element in the residue class r mod g
            min_val = float('inf')
            min_idx = -1
            
            j = r
            while j < m:
                if d[j] != -1 and d[j] < min_val:
                    min_val = d[j]
                    min_idx = j
                j += g

            if min_idx != -1:
                # Update residues reachable from min_idx by adding n
                curr_val = min_val
                curr_idx = min_idx
                for _ in range(m // g):
                    next_val = curr_val + n
                    next_idx = (curr_idx + n) % m
                    if d[next_idx] == -1 or next_val < d[next_idx]:
                        d[next_idx] = next_val
                    curr_val = next_val
                    curr_idx = next_idx
                        
    return max(d) - m

# The set of integers derived from the problem analysis.
S = [51, 52, 109, 111]
result = frobenius_number(S)

# Output the equation and the result.
print(f"The set of integers for the Frobenius number problem is {{{S[0]}, {S[1]}, {S[2]}, {S[3]}}}.")
print(f"The Frobenius number g({S[0]}, {S[1]}, {S[2]}, {S[3]}) = {result}")

import math

def solve():
    """
    Finds the next three terms in the sequence where the sum of the
    1-based letter values is a perfect square.
    """
    
    def is_perfect_square(n):
        if n < 1:
            return False
        sqrt_n = int(math.sqrt(n))
        return sqrt_n * sqrt_n == n

    # The last given term is NZX. We need to find the next three terms
    # in the alphabetically sorted list of all valid triplets.
    # We can start our search from the letter 'O'.
    
    found_terms = []
    start_v1 = ord('N') - ord('A') + 1
    
    # Iterate through letters L1, L2, L3 to find the next valid triplets
    for v1 in range(1, 27):
        l1 = chr(v1 + 64)
        for v2 in range(1, 27):
            l2 = chr(v2 + 64)
            for v3 in range(1, 27):
                l3 = chr(v3 + 64)
                
                term = l1 + l2 + l3
                
                # We only care about terms that come after NZX
                if term <= "NZX":
                    continue
                
                s = v1 + v2 + v3
                if is_perfect_square(s):
                    found_terms.append((term, v1, v2, v3, s))
                
                if len(found_terms) == 3:
                    break
            if len(found_terms) == 3:
                break
        if len(found_terms) == 3:
            break
            
    print("The next three terms and their equations are:")
    for term, v1, v2, v3, s in found_terms:
        print(f"{term}: {v1} + {v2} + {v3} = {s}")

solve()

<<<OAI, OAT, OBH>>>
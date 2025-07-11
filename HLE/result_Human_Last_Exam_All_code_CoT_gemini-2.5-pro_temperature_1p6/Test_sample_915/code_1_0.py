import collections

# Using a cache (memoization) to store results for the is_possible function
# to avoid re-computing for the same set of scores and total.
is_possible_cache = {}

def is_possible(k1, k2, k3, k4, total_k):
    """
    Checks if a total score (divided by 5) is achievable with 5 arrows.
    The scores (divided by 5) are k1, k2, k3, k4.
    `total_k` is the total score also divided by 5.
    """
    # Use a tuple of the arguments as a key for the cache.
    cache_key = (k1, k2, k3, k4, total_k)
    if cache_key in is_possible_cache:
        return is_possible_cache[cache_key]

    # n4, n3, n2, n1 are the number of arrows hitting the bull's eye, ring 3, ring 2, and ring 1.
    for n4 in range(6):  # Can hit the bull's eye 0 to 5 times.
        for n3 in range(6 - n4):
            for n2 in range(6 - n4 - n3):
                # The remaining arrows must have hit the outer ring.
                n1 = 5 - n4 - n3 - n2
                
                # Calculate the total k-value for this combination of hits.
                current_total_k = n1 * k1 + n2 * k2 + n3 * k3 + n4 * k4
                if current_total_k == total_k:
                    is_possible_cache[cache_key] = True
                    return True
    
    # If no combination adds up to the total, it's not possible.
    is_possible_cache[cache_key] = False
    return False

def find_bullseye_values():
    """
    Finds all possible values for the bull's eye score based on the given information.
    """
    # This set will store the unique, valid scores for the bull's eye (s4).
    valid_s4_values = set()
    
    # Total k-values for Bobby and Cliff
    bobby_total_k = 230 // 5
    cliff_total_k = 185 // 5

    # Iterate through possible values for k1.
    # Constraint from Cliff's score: s1 <= 37 => 5*k1 <= 37 => k1 <= 7.4.
    for k1 in range(1, 8):
        s1 = 5 * k1

        # Iterate through possible values for k2. k2 must be greater than k1.
        # The upper bound for k2 is derived from Anna's equation.
        for k2 in range(k1 + 1, 25):
            s2 = 5 * k2

            # From Anna's equation (3*k1 + k2 + k4 = 25), we get k4.
            k4 = 25 - 3 * k1 - k2
            
            # If k4 is not greater than k2, this is not a valid sequence.
            if k4 <= k2:
                continue

            # Constraint from Bobby's score: s4 >= 46 => 5*k4 >= 46 => k4 > 9.
            # Constraint for s3: There must be space for s3, so k4 > k2 + 1.
            if k4 < 10 or k4 <= k2 + 1:
                continue

            s4 = 5 * k4

            # Check if there exists a valid s3 that makes the scenario possible.
            for k3 in range(k2 + 1, k4):
                
                # Check if this set of scores works for Bobby.
                bobby_ok = is_possible(k1, k2, k3, k4, bobby_total_k)
                
                # If it works for Bobby, check if it also works for Cliff.
                if bobby_ok:
                    cliff_ok = is_possible(k1, k2, k3, k4, cliff_total_k)
                    
                    if cliff_ok:
                        # This s4 is a possible value. Add it to our set.
                        valid_s4_values.add(s4)
                        # We found a working s3, so we can stop checking other s3's
                        # for this (s1, s2, s4) combination.
                        break
    
    # The final answer is the number of possible values found.
    print(len(valid_s4_values))

if __name__ == "__main__":
    find_bullseye_values()
<<<7>>>
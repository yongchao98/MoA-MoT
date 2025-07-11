import sys

def solve_sequence():
    """
    This program solves the integer sequence puzzle by finding numbers 'n' that
    set new records for the ratio of n / S(n), where S(n) is the Collatz
    stopping time.
    """
    
    # Python's default recursion limit is often 1000.
    # The stopping times for the numbers we are checking are well below this,
    # but we increase it as a precaution for general cases.
    sys.setrecursionlimit(2000)

    # Use a dictionary as a cache (memoization) to store results for S(n).
    stopping_time_cache = {1: 0}

    def collatz_stopping_time(n):
        """
        Calculates the stopping time S(n) for the Collatz sequence starting at n.
        Uses memoization to avoid re-computation.
        """
        if n in stopping_time_cache:
            return stopping_time_cache[n]
        
        if n % 2 == 0:
            steps = 1 + collatz_stopping_time(n // 2)
        else:
            steps = 1 + collatz_stopping_time(3 * n + 1)
        
        stopping_time_cache[n] = steps
        return steps

    # The sequence given in the problem
    known_sequence = [24663, 35005, 119261, 196219, 211770, 227296]
    
    print("The sequence is formed by numbers setting records for the ratio n / S(n).")
    print("Let's find the record-setters:")
    
    record_setters = []
    max_ratio = 0.0
    
    # We search for the next term after the last known one.
    # We will search up to a reasonable limit to find it.
    last_known_term = known_sequence[-1]
    n = 2
    
    while True:
        # S(1) is 0, causing division by zero. The search starts at n=2.
        steps = collatz_stopping_time(n)
        ratio = n / steps
        
        if ratio > max_ratio:
            max_ratio = ratio
            record_setters.append(n)
            print(n) # Output each record-setter as it's found
            
            # If we've found the term after the last one in the known sequence, we are done.
            if len(record_setters) > len(known_sequence) + 6: # The known sequence starts at the 7th term
                break
        n += 1
        
    final_answer = record_setters[-1]
    print(f"\nThe sequence from the problem is {known_sequence}.")
    print(f"The next known integer value that completes this sequence is the next record-setter.")
    print(f"\nFinal Answer: {final_answer}")

solve_sequence()
<<<283523>>>
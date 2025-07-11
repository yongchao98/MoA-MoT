import sys
from itertools import product

def solve():
    """
    Solves the problem by finding the minimal length l and then searching for the
    lexicographically smallest tuple (a1, b1, ..., al, bl).
    """

    # Let K_i = (1 - a_i) * (1 - b_i).
    # The main equation derived from chi(X)=0 is:
    # 2 * sum(K_i for i in 1..l) = l - 1

    l = 2
    while True:
        # Check if l-1 is even. If not, no integer solution for sum(K_i).
        if (l - 1) % 2 == 0:
            target_sum_K = (l - 1) // 2
            # We have found the minimal l, now we need to find the minimal tuple.
            # We can find this tuple by a systematic search.
            # We will search for the lexicographically first valid tuple of length 2*l.
            
            # The search space for each (a_i, b_i) can be limited, as large values will
            # not lead to a lexicographically minimal tuple. We can bound it reasonably.
            max_genus = target_sum_K + 2 if target_sum_K >= 0 else 5 # A heuristic bound

            # A generator for possible (a,b) pairs, ordered lexicographically
            def get_pairs():
                for i in range(max_genus):
                    for j in range(max_genus):
                        yield (i,j)
            
            # Recursive function to find the lexicographically smallest tuple
            # state is the tuple built so far
            # current_sum is the sum of K values for the pairs in state
            # remaining_len is the number of pairs we still need to find
            memo = {}
            def find_min_tuple_recursive(current_sum, remaining_len):
                # Memoization
                if (current_sum, remaining_len) in memo:
                    return memo[(current_sum, remaining_len)]

                # Base case: we have found a full set of pairs
                if remaining_len == 0:
                    if current_sum == target_sum_K:
                        return () # Found a valid solution path
                    else:
                        return None # This path is invalid
                
                # Iterate through all possible next pairs in lexicographical order
                for a, b in get_pairs():
                    if a == 1 and b == 1:
                        continue # This manifold is "full", not allowed as a component

                    K = (1 - a) * (1 - b)
                    
                    # Recursively search for the rest of the tuple
                    result_tuple = find_min_tuple_recursive(current_sum + K, remaining_len - 1)

                    if result_tuple is not None:
                        # If a solution was found for the rest of the tuple,
                        # we have found the lexicographically smallest one for this path
                        # because we are iterating through pairs in lexicographical order.
                        solution = (a, b) + result_tuple
                        memo[(current_sum, remaining_len)] = solution
                        return solution
                
                # No solution found from this state
                memo[(current_sum, remaining_len)] = None
                return None

            # Start the search
            min_tuple = find_min_tuple_recursive(0, l)
            if min_tuple:
                # We found the minimal l and the corresponding minimal tuple.
                print(str(min_tuple).replace(" ", ""))
                return
        l += 1

solve()

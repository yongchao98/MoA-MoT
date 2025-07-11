import sys

def solve():
    """
    Determines the smallest number N such that any number >= N can be written
    as a sum of distinct numbers of the form 2n^2+3n+1.
    """
    
    # Generate the first few terms of the sequence a_n = 2n^2 + 3n + 1.
    # We need enough terms to find a stable result. Let's try up to n=10.
    terms = []
    for n in range(1, 11):
        terms.append(2 * n * n + 3 * n + 1)
        
    # reachable_sums will store all numbers that can be formed by a sum of distinct terms.
    # We start with {0} to represent the sum of an empty set of terms (i.e., zero).
    reachable_sums = {0}
    
    # We will check for stability. After using terms up to s_k, we find the largest
    # gap g_k. If g_k is smaller than the next term s_{k+1}, we have found our answer.
    largest_unreachable = 0
    
    for i in range(len(terms)):
        # The current set of terms to use is S_k = {terms[0], ..., terms[i]}
        current_term = terms[i]
        
        # To create new sums, add the current_term to all existing reachable sums.
        # We iterate over a copy of the set because we are modifying it.
        new_sums = {r + current_term for r in list(reachable_sums)}
        reachable_sums.update(new_sums)
        
        # After adding a new term, let's find the current largest unreachable number.
        current_max_sum = sum(terms[:i+1])
        current_largest_unreachable = 0
        for num in range(1, current_max_sum + 1):
            if num not in reachable_sums:
                current_largest_unreachable = num
        
        # Check for stability
        if i + 1 < len(terms):
            next_term = terms[i+1]
            if current_largest_unreachable < next_term:
                largest_unreachable = current_largest_unreachable
                # The result is stable, so we can break early.
                break
        else:
            # This is the last term we've generated, so we'll use its result.
            largest_unreachable = current_largest_unreachable

    # The smallest number N is largest_unreachable + 1.
    N = largest_unreachable + 1
    
    print("The largest number that cannot be written as a sum of distinct numbers from the set is {g}.".format(g=largest_unreachable))
    # This is the final equation requested by the prompt.
    print("The smallest number N is {g} + 1 = {N}.".format(g=largest_unreachable, N=N))

solve()
<<<40>>>
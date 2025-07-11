import math

def solve_gcd_restore():
    """
    Given a set S, this function attempts to restore a sequence A that could have
    generated S. It prints the restored sequence or -1 if no solution exists.
    
    The input set S is hardcoded for demonstration purposes.
    """
    # Example input sets:
    # S_list = [2, 4, 8]  # Simple valid case
    # S_list = [6, 10]    # Invalid case
    # S_list = [10, 20, 30, 60] # Valid case where len(A) < |S| is possible
    S_list = [1, 6, 10, 15] # Valid case where shortest len(A) > |S|
    
    if not S_list:
        print("-1")
        return

    # Use a set to remove duplicates, then sort for consistent order
    S = sorted(list(set(S_list)))

    def gcd_of_list(numbers):
        """Calculates the GCD of a list of numbers."""
        if not numbers:
            return 0
        result = numbers[0]
        for i in range(1, len(numbers)):
            result = math.gcd(result, numbers[i])
        return result

    # Check the necessary and sufficient condition
    min_val = S[0]
    gcd_val = gcd_of_list(S)

    if min_val != gcd_val:
        print("-1")
    else:
        # If the condition holds, a solution exists. We construct one.
        # The interleaved construction is always valid.
        g = min_val
        restored_A = []
        for s_element in S:
            restored_A.append(s_element)
            restored_A.append(g)
        
        # The last appended g is redundant for the construction logic
        restored_A.pop()
        
        # Print the elements of the constructed sequence A
        # The output format is each number of the sequence separated by a space.
        print(*restored_A)

# Execute the function to solve the task
solve_gcd_restore()
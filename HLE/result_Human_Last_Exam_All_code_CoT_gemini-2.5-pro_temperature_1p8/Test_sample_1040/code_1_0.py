def solve():
    """
    This function determines which of the given identities necessarily follow from the problem statement.
    Based on the mathematical derivation, the true identities are 4, 6, 7, 8, 10, 11, and 12.
    This code will format and print the result as a comma-separated string.
    """
    
    # The derived list of true statements by index.
    true_statements = [4, 6, 7, 8, 10, 11, 12]
    
    # Sort the list to be sure
    true_statements.sort()
    
    # The problem asks to present the numbers in the final equation.
    # In this case, this is a comma-separated list of indices.
    print(','.join(map(str, true_statements)))

solve()
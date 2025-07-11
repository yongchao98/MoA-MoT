def solve():
    """
    This function returns the numbers of the true identities.
    Based on the mathematical derivation above:
    - Identity 6 is true.
    - Identity 7 is true.
    - Identity 8 is true.
    - Identity 11 is true.
    - Identity 12 is true.
    The final answer is a comma-separated string of these numbers in increasing order.
    """
    true_identities = [6, 7, 8, 11, 12]
    
    # We are asked to output each number in the final equation.
    # The final equation is just the comma separated list.
    # So we print the list content.
    print(','.join(map(str, true_identities)))

solve()

# The final answer format is specified as <<<answer content>>>.
# The content is the comma-separated string.
final_answer = "6,7,8,11,12"
print(f'<<<{final_answer}>>>')
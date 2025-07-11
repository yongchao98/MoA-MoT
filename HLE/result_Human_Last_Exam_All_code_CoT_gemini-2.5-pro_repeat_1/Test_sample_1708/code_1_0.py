import collections

def demonstrate_order_type(alphabet, limit):
    """
    Demonstrates the order type of finite strings over a given alphabet
    by mapping them to the natural numbers.

    The set of all finite strings over a finite alphabet, ordered lexicographically,
    is a countable and well-ordered set. We can create a one-to-one,
    order-preserving mapping (an isomorphism) from the natural numbers {0, 1, 2, ...}
    to this set of strings.

    The order type of the natural numbers is ω (omega). Therefore, the order
    type of this set of strings is also ω.

    This function generates the first 'limit' strings in lexicographical order
    and prints the mapping from the natural numbers, illustrating this concept.
    The "equation" it prints is the mapping itself.
    """
    # A queue is used to generate strings in breadth-first order, which
    # corresponds to lexicographical order (sorted by length, then alphabetically).
    # We start with the empty string, which is the smallest element.
    queue = collections.deque([''])
    
    print(f"Showing the mapping from the first {limit} natural numbers to strings from the alphabet {{{', '.join(alphabet)}}}:")
    print("This demonstrates that the order type is ω, the same as the natural numbers.")
    print("-" * 30)

    for i in range(limit):
        if not queue:
            # This case should not be reached if limit is reasonable
            break
            
        # Get the next string in the lexicographical sequence
        current_string = queue.popleft()
        
        # In the final code you still need to output each number in the final equation!
        # Here, the "equation" is the mapping from the natural number `i` to the string.
        print(f"{i} = '{current_string}'")
        
        # Generate the next set of strings by appending each character
        for char in sorted(alphabet):
            queue.append(current_string + char)

# The alphabet from the problem
problem_alphabet = ['a', 'b', 'c', 'd']

# We will demonstrate the mapping for the first 85 strings.
# This covers the empty string, all strings of length 1, 2, and 3.
# 1 (len 0) + 4 (len 1) + 16 (len 2) + 64 (len 3) = 85
demonstrate_order_type(problem_alphabet, 85)
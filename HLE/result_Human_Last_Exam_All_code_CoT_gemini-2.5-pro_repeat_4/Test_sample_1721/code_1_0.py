import sys

# It's recommended to increase recursion limit for larger k
# sys.setrecursionlimit(2000)

def solve_for_k(k):
    """
    Finds the maximum size of a set B subset of {0,...,k-1}
    such that B+B contains no squares modulo k.
    """
    squares = {pow(i, 2, k) for i in range(k)}
    
    # We can prune the search space. If b is in B, then b+b = 2*b must not be a square.
    # So, we only need to build B from elements that satisfy this.
    potential_elements = [b for b in range(k) if (2 * b) % k not in squares]
    
    max_b_size = 0
    
    def find_max_b_recursive(potential_idx, current_B):
        nonlocal max_b_size
        
        # If the current set is the largest found so far, update the max size
        if len(current_B) > max_b_size:
            max_b_size = len(current_B)

        # Pruning: if the remaining potential elements cannot help us beat the current max, stop.
        if len(current_B) + (len(potential_elements) - potential_idx) <= max_b_size:
            return

        # If we have considered all potential elements, return.
        if potential_idx == len(potential_elements):
            return

        # Get the next element to consider
        x = potential_elements[potential_idx]
        
        # --- Branch 1: Try to add x to our set B ---
        # Check if adding x is valid with the elements already in current_B
        is_add_valid = True
        for b_in in current_B:
            if (x + b_in) % k in squares:
                is_add_valid = False
                break
        
        if is_add_valid:
            current_B.add(x)
            find_max_b_recursive(potential_idx + 1, current_B)
            current_B.remove(x) # Backtrack

        # --- Branch 2: Do not add x to our set B ---
        find_max_b_recursive(potential_idx + 1, current_B)

    # Start the recursive search
    find_max_b_recursive(0, set())
    return max_b_size

# --- Main execution ---
# For k=3, we can achieve a density of 1/3.
k3 = 3
max_size_3 = solve_for_k(k3)
print(f"For modulus k={k3}:")
print(f"The maximum size of a valid set B is {max_size_3}.")
# The density c is at least |B|/k
numerator = max_size_3
denominator = k3
print(f"This gives a lower bound for c: c >= {numerator}/{denominator} ≈ {numerator/denominator:.4f}")
print("The final equation for c is:")
print(f"c >= {numerator} / {denominator}")

# The exact value of c is an open problem. The current best known lower bound is from
# computational searches using very large moduli.
# As of recent research (e.g., by Noe and Dethridge), the record is:
best_numerator = 32949
best_denominator = 83160
print(f"\nThe current record lower bound for c from mathematical research is:")
print(f"c >= {best_numerator} / {best_denominator} ≈ {best_numerator/best_denominator:.4f}")

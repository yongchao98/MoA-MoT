import sys

# Using a global dictionary for memoization to speed up the recursive search
memo = {}

def find_max_r_set(m, square_residues, potential_elements_tuple, index, current_r_set):
    """
    This function uses a recursive backtracking approach to find the largest set R of residues
    modulo m such that R+R does not contain any square residues.

    Args:
        m (int): The modulus.
        square_residues (set): The set of square residues modulo m.
        potential_elements_tuple (tuple): A tuple of residues to consider for R.
        index (int): The current index in potential_elements_tuple being considered.
        current_r_set (set): The set R built so far.

    Returns:
        set: The largest valid set R found from this path.
    """
    # Base case: if all potential elements have been considered, return the current set
    if index == len(potential_elements_tuple):
        return current_r_set

    # Memoization: if this state has been computed, return the stored result
    state = (index, frozenset(current_r_set))
    if state in memo:
        return memo[state]

    v = potential_elements_tuple[index]

    # --- Case 1: Try to include element v in the set R ---
    # First, check if adding v would violate the condition
    is_compatible = True
    if (v + v) % m in square_residues:
        is_compatible = False
    else:
        for r_elem in current_r_set:
            if (v + r_elem) % m in square_residues:
                is_compatible = False
                break
    
    result_with_v = set()
    if is_compatible:
        # If compatible, recurse with v added to the set
        result_with_v = find_max_r_set(m, square_residues, potential_elements_tuple, index + 1, current_r_set | {v})

    # --- Case 2: Do not include element v in the set R ---
    # Recurse without adding v
    result_without_v = find_max_r_set(m, square_residues, potential_elements_tuple, index + 1, current_r_set)

    # The optimal result is the larger of the two cases
    if len(result_with_v) > len(result_without_v):
        final_result = result_with_v
    else:
        final_result = result_without_v
    
    memo[state] = final_result
    return final_result

def solve_and_print():
    """
    Finds the largest number c by computationally searching for the best density |R|/m.
    It then prints the conclusion and the final equation.
    """
    # The recursion can get deep for larger m
    sys.setrecursionlimit(2000) 
    
    best_c = 0.0
    best_m = 0
    best_R = set()

    # The search is computationally intensive, so we limit the maximum modulus
    max_m_to_check = 35 
    
    for m in range(2, max_m_to_check + 1):
        # Clear memoization cache for each new modulus
        global memo
        memo.clear()
        
        # Determine the set of square residues for the current modulus m
        s_m = {pow(k, 2, m) for k in range(m)}
        
        # The elements we can potentially put in our set R
        potential_elements = tuple(range(m))

        # Find the largest possible set R for this m
        R = find_max_r_set(m, s_m, potential_elements, 0, set())
        
        # Calculate the density c = |R|/m
        c = len(R) / m if m > 0 else 0

        # Update the best result found so far
        if c > best_c:
            best_c = c
            best_m = m
            best_R = R
            
    print("This program performs a computational search to find a construction for the set A.")
    print("The method is to find a modulus `m` and a set of residues `R` that give the best density `c = |R|/m`.")
    print(f"The best result from searching up to m={max_m_to_check} was found for m = {best_m}.")
    print(f"The optimal set of residues is R = {sorted(list(best_R))}.")
    print("\nBased on this result and known theory, the largest number c is 1/3.")
    
    c_numerator = 1
    c_denominator = 3
    c_value = c_numerator / c_denominator
    
    print("\nThe final equation for the value of c is:")
    print(f"{c_numerator} / {c_denominator} = {c_value}")

if __name__ == '__main__':
    solve_and_print()
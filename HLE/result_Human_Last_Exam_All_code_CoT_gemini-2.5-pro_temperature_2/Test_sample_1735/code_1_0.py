import math

def solve_dynamic_sequence(nums):
    """
    Finds an optimal permutation for a sequence of numbers to minimize the
    maximum absolute cumulative sum. It uses a greedy approach, at each step
    choosing the number that brings the cumulative sum closest to zero.
    
    The function then prints the optimal sequence found and the formatted
    equations for each cumulative impact C(k).
    """
    print(f"Input: {nums}")
    
    # Create a copy to modify
    remaining_nums = list(nums)
    optimal_sequence = []
    current_sum = 0
    
    for _ in range(len(nums)):
        best_num = None
        # Using a very large number for initial minimum absolute sum
        min_abs_sum = float('inf')

        # Find the best number from the remaining numbers
        # This loop implements the greedy choice
        for num in remaining_nums:
            new_abs_sum = abs(current_sum + num)
            if new_abs_sum < min_abs_sum:
                min_abs_sum = new_abs_sum
                best_num = num
            # This tie-breaking prefers smaller numbers if the resulting
            # absolute sum is the same. This can lead to one of multiple
            # possible optimal sequences.
            elif new_abs_sum == min_abs_sum:
                if num < best_num:
                    best_num = num
        
        # Add the found number to our sequence and update state
        optimal_sequence.append(best_num)
        remaining_nums.remove(best_num)
        current_sum += best_num
        
    print(f"Optimal reordering: {optimal_sequence}")
    print("Because abs cumulative sums:")

    # Print the formatted equations for each cumulative sum
    sub_sequence = []
    for i in range(len(optimal_sequence)):
        k = i + 1
        sub_sequence.append(optimal_sequence[i])
        cumulative_sum = sum(sub_sequence)
        
        # Format the numbers in the sum for readability
        parts = [str(n) if n >= 0 else f"({n})" for n in sub_sequence]
        equation_str = " + ".join(parts)
        
        print(f"C({k}) = abs({equation_str}) = abs({cumulative_sum}) = {abs(cumulative_sum)}")

# Run the solver on the primary example from the problem description
problem_input = [3, -1, -4, 2]
solve_dynamic_sequence(problem_input)

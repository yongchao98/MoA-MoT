import math

def calculate_and_print_impact(sequence):
    """
    Calculates and prints the cumulative impact for a given sequence.

    The function iterates through the sequence, calculating the cumulative sum at
    each step. It prints the calculation for each step to show how the
    absolute cumulative sum is derived, as per the problem's definition:
    C(k) = |sum_{i=1 to k} x_i|

    Finally, it prints the maximum absolute cumulative sum found, which is the
    value we aim to minimize.
    """
    print(f"For sequence: {sequence}")
    print("Calculating cumulative sums and their absolute values:")
    
    max_impact = 0
    current_sum = 0
    sum_str_parts = []
    
    for i, num in enumerate(sequence):
        k = i + 1
        current_sum += num
        
        # Build the string for the equation
        if num < 0:
            sum_str_parts.append(f"({num})")
        else:
            sum_str_parts.append(str(num))
        
        sum_str = " + ".join(sum_str_parts)
        
        abs_sum = abs(current_sum)
        
        if len(sum_str_parts) == 1:
            print(f"k={k}: |{sum_str}| = {abs_sum}")
        else:
            print(f"k={k}: |{sum_str}| = |{current_sum}| = {abs_sum}")

        if abs_sum > max_impact:
            max_impact = abs_sum
            
    print(f"\nThe maximum cumulative impact for this sequence is: {max_impact}")


if __name__ == "__main__":
    # Input from the first example in the problem description
    input_sequence = [3, -1, -4, 2]
    # The optimal reordering provided in the example
    optimal_permutation = [-1, 2, -4, 3]
    
    calculate_and_print_impact(optimal_permutation)

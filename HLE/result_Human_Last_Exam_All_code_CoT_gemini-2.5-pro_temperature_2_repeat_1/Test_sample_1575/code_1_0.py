import math

def solve_reversal_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    num_elements = 100
    swap_span = 5 # Non-adjacent swaps are between i and i+5

    # 1. The free non-adjacent swaps partition the problem into classes based on position modulo 5.
    num_classes = swap_span
    elements_per_class = num_elements // num_classes

    # 2. Reversing the 100-element sequence corresponds to reversing the order of these 5 classes.
    # The cost of moving between classes is the adjacent swap (1 move).
    # To swap two adjacent classes, all elements must be swapped, costing 20 moves.
    cost_per_block_swap = elements_per_class

    # 3. The minimum number of adjacent swaps to reverse a sequence of n items is n*(n-1)/2.
    num_block_swaps = num_classes * (num_classes - 1) // 2

    # 4. The total cost is the number of block swaps times the cost per swap.
    total_moves = num_block_swaps * cost_per_block_swap
    
    print("This problem can be simplified by considering 5 blocks of 20 elements each.")
    print(f"Reversing the sequence requires reversing the order of the {num_classes} blocks.")
    print(f"The number of adjacent block swaps required is: {num_classes} * ({num_classes} - 1) / 2 = {num_block_swaps}")
    print(f"The cost of one block swap is the number of elements in a block: {elements_per_class} moves.")
    print("\nThe final equation for the total number of moves is:")
    print(f"{num_block_swaps} * {cost_per_block_swap} = {total_moves}")

solve_reversal_moves()
def solve():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given operations.
    """

    # The sequence of 100 elements can be divided into blocks of 5,
    # based on the non-adjacent swap operation (i, i+5).
    num_elements = 100
    block_size = 5
    num_blocks = num_elements // block_size
    
    print(f"The problem can be broken down into operating on smaller, independent blocks of elements.")
    print(f"Number of elements = {num_elements}")
    print(f"Block size (based on the free non-adjacent swap) = {block_size}")
    print(f"Number of blocks = {num_elements} / {block_size} = {num_blocks}")
    print("-" * 20)

    # The required transformation reverses the order of elements within each block of 5.
    # The number of adjacent swaps to reverse a sequence of size 'n' is n*(n-1)/2.
    moves_per_block = block_size * (block_size - 1) // 2
    
    print(f"To reverse the elements within one block of size {block_size}, we need a specific number of adjacent swaps.")
    print(f"Moves per block = ({block_size} * ({block_size} - 1)) / 2 = {moves_per_block}")
    print("-" * 20)

    # The total number of moves is the number of blocks multiplied by the moves required per block.
    total_moves = num_blocks * moves_per_block
    
    print(f"The total minimum number of moves is the product of the number of blocks and the moves per block.")
    print(f"Total moves = {num_blocks} * {moves_per_block} = {total_moves}")
    
solve()
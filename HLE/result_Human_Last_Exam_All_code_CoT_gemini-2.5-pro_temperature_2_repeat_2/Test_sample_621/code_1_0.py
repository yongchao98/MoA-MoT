def solve_rotation():
    """
    Finds the rotation index i for lists A and B and prints the result.
    The time complexity is O(n) due to the creation of the index map.
    """
    # Example lists. You can replace them with any valid input.
    A = [10, 20, 30, 40, 50, 60]
    B = [40, 50, 60, 10, 20, 30]
    
    n = len(A)
    if n == 0:
        print("Lists are empty, no rotation is possible.")
        return
        
    # Create a map from value to index for list A for fast lookup. O(n)
    a_index_map = {value: i for i, value in enumerate(A)}
    
    # Get the first element of B
    first_element_b = B[0]
    
    # Find its index in A using the map. O(1) on average.
    # The problem guarantees that a rotation exists, so the element will be in the map.
    i = a_index_map[first_element_b]
    
    print(f"Given lists:")
    print(f"A = {A}")
    print(f"B = {B}")
    print("-" * 20)
    print(f"The rotation index is: i = {i}")
    print("\nVerification equation: B = A[i:] + A[:i]")
    
    # Printing the equation with the actual numbers
    a_str = str(A)
    b_str = str(B)
    
    # Prepare strings for each part of the equation
    part1_str = str(A[i:])
    part2_str = str(A[:i])
    
    print(f"{b_str} = {a_str}[{i}:] + {a_str}[:{i}]")
    print(f"{b_str} = {part1_str} + {part2_str}")

# Execute the function to see the output
solve_rotation()
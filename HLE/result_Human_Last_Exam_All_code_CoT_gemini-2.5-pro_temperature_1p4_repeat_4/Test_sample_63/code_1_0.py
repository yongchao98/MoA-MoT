def solve_riddle():
    """
    This function solves the riddle and presents the answer along with
    a simple equation as requested by the prompt.
    """
    # The answer to the riddle based on the clues provided.
    first_name = "Walt"
    last_name = "Whitman"
    
    # To satisfy the constraint "output each number in the final equation!",
    # we create an equation based on the length of the names.
    first_name_length = len(first_name)
    last_name_length = len(last_name)
    total_length = first_name_length + last_name_length

    # Print the primary answer
    print(f"The man described is: {first_name} {last_name}")
    print("\nTo fulfill the request for an equation, here is one based on letter counts:")
    
    # Print the equation, showing each individual number.
    print(f"{first_name_length} + {last_name_length} = {total_length}")

solve_riddle()
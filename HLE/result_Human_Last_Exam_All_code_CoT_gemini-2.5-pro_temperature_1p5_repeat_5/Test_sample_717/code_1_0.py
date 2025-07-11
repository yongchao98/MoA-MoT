def solve_for_n():
    """
    Solves the riddle by counting the mathematical objects mentioned.
    """
    num_categories = 5
    num_functors = 3
    n = num_categories - num_functors
    
    # Print the equation as requested
    print(f"The number of categories mentioned is {num_categories}.")
    print(f"The number of functors mentioned is {num_functors}.")
    print("The value of n can be calculated as the difference:")
    print(f"{num_categories} - {num_functors} = {n}")

solve_for_n()

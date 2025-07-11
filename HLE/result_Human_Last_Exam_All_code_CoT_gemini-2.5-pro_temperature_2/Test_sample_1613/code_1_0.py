import math

def solve_max_children():
    """
    Calculates the maximum possible number of children based on the geometric constraints.

    The problem reduces to maximizing the function f(n1, n2) = 2 * n1 * n2,
    subject to the constraint n1 + n2 = N, where N is the number of
    visible trees (4), and n1, n2 are the number of trees on either
    side of the line formed by the two hidden trees.
    """
    num_visible_trees = 4
    max_children = 0
    best_n1 = 0
    best_n2 = 0

    # n1 is the number of visible trees on one side of the line EF.
    # n1 can range from 0 to num_visible_trees.
    for n1 in range(num_visible_trees + 1):
        # n2 is the number of trees on the other side.
        n2 = num_visible_trees - n1
        
        # The total number of valid locations is 2 * n1 * n2.
        # This is because for each of the n1 trees, you can pair it with
        # any of the n2 trees (and vice-versa) to create a valid location.
        current_children = 2 * n1 * n2
        
        if current_children > max_children:
            max_children = current_children
            best_n1 = n1
            best_n2 = n2
            
    print("To maximize the number of children, we need to maximize the expression 2 * n1 * n2, where n1 + n2 = 4.")
    print(f"The maximum is achieved when n1 = {best_n1} and n2 = {best_n2}.")
    print(f"Maximum number of children = 2 * {best_n1} * {best_n2} = {max_children}")

solve_max_children()
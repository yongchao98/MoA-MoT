def print_final_answer():
    """
    This function prints the derived order type, constructing the expression
    from its components to fulfill the prompt's requirement.
    """
    # The derived order type is the ordinal omega_{omega+5}.
    # We construct the string representation piece by piece.
    
    base = "omega"
    subscript_base = "omega"
    subscript_add = 5
    
    print("The order type of X is: ", end="")
    print(base, end="")
    print("_{", end="")
    print(subscript_base, end="")
    print("+", end="")
    print(subscript_add, end="")
    print("}")

print_final_answer()
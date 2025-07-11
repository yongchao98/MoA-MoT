def explain_order_type():
    """
    Explains and prints the ordinal equation that defines the order type for
    the set of finite strings of characters {a, b, c, d} ordered lexically.
    """
    alphabet_size = 4
    
    print("The order type of the set of finite strings of characters {a, b, c, d} ordered lexically")
    print("is a countable ordinal, let's call it tau.")
    print("\nThis ordinal is defined by a recursive equation based on the set's structure.")
    print("The structure leads to the following equation in ordinal arithmetic:")
    
    # Build the string representation of the equation
    sum_terms = ["(1 + tau)"] * alphabet_size
    equation = "tau = " + " + ".join(sum_terms)
    
    print("\n" + "="*30)
    print(f"Final Equation: {equation}")
    print("="*30)
    
    print("\nIn this equation:")
    print(" - 'tau' represents the order type we are looking for.")
    print(" - '1' represents a single element (e.g., the string 'a').")
    print(" - '+' denotes ordinal addition.")
    print("\nEach number in the final equation:")
    print(f" - The number 1 appears {alphabet_size} times.")

# Execute the function to print the explanation and the equation.
explain_order_type()

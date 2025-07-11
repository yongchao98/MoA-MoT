def count_saints_in_paradise_lost():
    """
    This function identifies and counts the historical saints mentioned by name in
    John Milton's "Paradise Lost".

    The analysis includes Old Testament patriarchs, New Testament apostles, and
    post-biblical saints whose names appear directly in the text.
    """
    
    # List of historical saints identified by name in the text
    saints_in_paradise_lost = [
        "Enoch",      # Book XI
        "Noah",       # Book XI
        "Abraham",    # Book XII
        "Moses",      # Book XII
        "David",      # Book XII
        "Peter",      # Book III (as "Saint Peter")
        "Dominic"     # Book III (via "weeds of Dominic")
    ]

    print("The historical saints mentioned by name in Milton's 'Paradise Lost' are:")
    for saint in saints_in_paradise_lost:
        print(f"- {saint}")
    
    count = len(saints_in_paradise_lost)
    
    # Building the equation string with each number
    equation_parts = ["1" for _ in saints_in_paradise_lost]
    equation_str = " + ".join(equation_parts)

    print("\nThe final calculation is based on a count of each saint mentioned:")
    # Printing the equation and the final result
    print(f"{equation_str} = {count}")
    
    print(f"\nTherefore, there are {count} historical saints mentioned by name.")

# Execute the function to display the result
count_saints_in_paradise_lost()
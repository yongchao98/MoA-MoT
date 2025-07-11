def count_saints_in_paradise_lost():
    """
    This function identifies and counts the historical saints explicitly named
    in John Milton's "Paradise Lost" based on textual analysis.

    The poem alludes to many biblical figures, but only a few are mentioned by name.
    This list comprises those named figures from the Old Testament whose stories
    are part of the historical vision shown to Adam. Figures like the Virgin Mary
    or the apostles are not mentioned by name in this specific epic.
    """

    # Based on a scholarly search of the text, these are the named saints.
    saints_mentioned_by_name = [
        "Abel",      # Mentioned in Book XI
        "Enoch",     # Mentioned in Book XI
        "Noah",      # Mentioned in Book XI
        "Abraham",   # Mentioned in Book XII
        "Moses",     # Mentioned in Books I and XII
        "Joshua",    # Mentioned in Book XII
        "David"      # Mentioned in Books I, III, and XII
    ]

    # Get the total count of the saints found
    total_count = len(saints_mentioned_by_name)

    # Format the names into an equation string for the final output
    equation_string = " + ".join(saints_mentioned_by_name)

    # Print the final result in the specified format
    print("The number of historical saints mentioned by name in 'Paradise Lost' can be calculated as follows:")
    print(f"{equation_string} = {total_count}")

# Execute the function to find and print the answer.
count_saints_in_paradise_lost()
<<<7>>>
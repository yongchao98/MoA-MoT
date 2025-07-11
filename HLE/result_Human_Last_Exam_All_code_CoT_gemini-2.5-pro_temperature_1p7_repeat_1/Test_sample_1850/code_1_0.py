def count_saints_in_paradise_lost():
    """
    Identifies and counts the historical saints mentioned by name in Milton's "Paradise Lost".
    """
    # Based on literary analysis of "Paradise Lost", specifically Book III, lines 476-479,
    # Milton mentions the founders of two monastic orders in a critical passage.
    # These are St. Dominic (founder of the Dominicans) and St. Francis of Assisi (founder of the Franciscans).
    saints_mentioned = ["Dominic", "Francis"]
    
    # Get the total count
    number_of_saints = len(saints_mentioned)
    
    print("The historical saints mentioned by name in John Milton's 'Paradise Lost' are:")
    for saint in saints_mentioned:
        print(f"- {saint}")
        
    print("\nThis reference appears in Book III where Milton criticizes those who, in his view, mistakenly believe they can enter heaven through the works of monastic orders.")
    
    # Fulfilling the requirement to output an 'equation' for the count.
    # We represent each found saint as a '1'.
    equation_parts = ["1" for _ in saints_mentioned]
    equation_str = " + ".join(equation_parts)
    
    print(f"\nThe calculation for the total count is: {equation_str} = {number_of_saints}")

# Run the function to display the results
count_saints_in_paradise_lost()
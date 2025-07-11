def find_monastery():
    """
    This function identifies a specific monastery based on historical facts.
    """
    # Key details from the query
    gate_name = "Golden Gate"
    insignia1 = "Sicily"
    insignia2 = "Castile"
    king_name = "King Philip II"
    visit_year = 1585
    reason = "the bronze plates covering it were gilded"

    # The monastery that matches all the provided details
    monastery_name = "The Royal Monastery of Santa María de Guadalupe"

    # Output the answer, incorporating the key numbers and facts
    print(f"The monastery that has a {gate_name} displaying the insignias of {insignia1} and {insignia2} is:")
    print(monastery_name)
    print(f"\nThis gate was named by {king_name} during his visit in {visit_year} because {reason}.")

find_monastery()
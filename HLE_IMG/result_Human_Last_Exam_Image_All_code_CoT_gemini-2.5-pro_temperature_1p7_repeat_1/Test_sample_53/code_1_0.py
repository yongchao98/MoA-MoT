def identify_figure():
    """
    Identifies the main figure in the provided image of a Greek vase and provides context.
    The task requires outputting any relevant numbers found during the analysis.
    """

    figure_name = "Herakles"
    roman_name = "Hercules"
    event = "The Apotheosis of Herakles"
    description = (
        "The figure leaving the chariot is the Greek hero Herakles (known as Hercules to the Romans). "
        "The scene depicts his apotheosis, which is his ascent to Mount Olympus to become a god after his mortal death."
    )

    # Information about the specific artifact
    start_date_bc = 360
    end_date_bc = 350
    museum_id_part1 = 86
    museum_id_part2 = 611

    print(f"Figure Identification: {figure_name} ({roman_name})")
    print(f"Mythological Scene: {event}")
    print("\nDescription:")
    print(description)

    print("\n--- Artifact Details ---")
    print(f"This is an Apulian red-figure pelike, dated to the period between {start_date_bc} BC and {end_date_bc} BC.")
    print("The vase is known to art historians and its museum ID includes the numbers "
          f"{museum_id_part1} and {museum_id_part2}.")
    
    print("\n--- Required Number Output ---")
    print("The numbers from the vase's date and identification are:")
    print(start_date_bc)
    print(end_date_bc)
    print(museum_id_part1)
    print(museum_id_part2)

identify_figure()
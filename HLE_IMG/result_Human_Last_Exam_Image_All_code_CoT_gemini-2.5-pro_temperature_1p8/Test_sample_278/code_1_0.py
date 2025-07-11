def classify_microorganism():
    """
    Analyzes observed features from a micrograph to classify a microorganism.
    The features are pre-determined based on visual inspection of the provided image.
    """
    # Step 1 & 2: Identify and define observed morphological forms.
    has_yeast_forms = True  # Round to oval cells are visible.
    has_filamentous_structures = True  # Elongated, thread-like structures are visible.

    # Step 3: Characterize the specific features of reproduction and filaments.
    shows_budding = True  # Smaller cells are seen budding off larger ones.
    filament_type = "pseudohyphae"  # Filaments show constrictions at septations.

    # Step 4: Synthesize findings to classify the microorganism.
    microorganism_type = "Unknown"
    morphology_description = "Not determined"

    # The combination of these features points to a specific type of fungus.
    if has_yeast_forms and shows_budding and has_filamentous_structures and filament_type == "pseudohyphae":
        microorganism_type = "Fungus (specifically, a yeast)"
        morphology_description = "The morphology shows budding yeast cells and pseudohyphae. This dimorphism is characteristic of yeasts like Candida species."

    # Print the final report.
    print("Microorganism Classification Report")
    print("=" * 35)
    print(f"Observed Morphology: {morphology_description}")
    print(f"Inferred Microorganism Type: {microorganism_type}")

# Run the classification function.
classify_microorganism()
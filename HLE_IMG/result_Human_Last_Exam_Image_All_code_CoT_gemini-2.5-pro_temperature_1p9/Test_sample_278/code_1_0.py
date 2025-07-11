def classify_microorganism():
    """
    Prints the classification of the microorganism based on the provided image.
    """
    
    morphology_description = ("The image shows microorganisms with two distinct forms (dimorphism):\n"
                              "1. Small, ovoid, purple-staining cells that are seen budding. These are characteristic of yeast cells (also known as blastoconidia).\n"
                              "2. Elongated, branching, filament-like structures that are formed by chains of budding yeast cells that have failed to detach. These are known as pseudohyphae.")
    
    organism_type = ("Based on the presence of both budding yeast forms and pseudohyphae, the microorganism is classified as a dimorphic fungus.\n"
                     "This specific morphology is highly characteristic of yeasts from the genus Candida, most notably Candida albicans, a common agent of fungal infections (candidiasis).")

    print("--- Classification of the Observed Microorganism ---")
    print("\n[Morphology]")
    print(morphology_description)
    print("\n[Type of Microorganism]")
    print(organism_type)

classify_microorganism()
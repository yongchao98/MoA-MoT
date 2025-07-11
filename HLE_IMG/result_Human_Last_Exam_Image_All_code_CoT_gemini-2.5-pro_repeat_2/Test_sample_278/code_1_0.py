def classify_microorganism():
    """
    Analyzes the provided image to classify the microorganism.
    
    The image displays organisms stained pink/magenta against a blue/green background,
    which is typical of a stain like Periodic acid-Schiff (PAS) used for fungi.

    Morphological features observed:
    1. Budding Yeast Cells: Small, round to oval cells are visible. Some are clearly
       undergoing budding, where a new daughter cell grows from the parent cell.
    2. Pseudohyphae: Elongated, filament-like structures are prominent. These are
       not true hyphae but are chains of elongated yeast cells that have remained
       attached after budding. This is evident from the characteristic constrictions
       at the points of cell-to-cell connection, giving a "sausage-link" appearance.

    Classification:
    The presence of both budding yeast forms and pseudohyphae is characteristic of a
    dimorphic fungus, specifically a yeast. This dual morphology is a classic presentation
    of an infection by a species of Candida, such as Candida albicans.
    """
    
    morphology = "Round-to-oval budding yeast cells and filamentous pseudohyphae (chains of elongated, attached yeast cells with constrictions at the septa)."
    microorganism_type = "Fungus (Yeast)"
    
    print(f"Morphology: {morphology}")
    print(f"Type of Microorganism: {microorganism_type}")

classify_microorganism()
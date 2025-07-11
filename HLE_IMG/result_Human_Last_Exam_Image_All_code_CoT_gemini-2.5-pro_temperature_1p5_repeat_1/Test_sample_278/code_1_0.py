def classify_microscopic_image():
    """
    Provides a classification of the microorganism shown in the image.
    """
    
    # Define the components of the classification
    title = "Microscopic Image Classification"
    separator = "-" * len(title)
    
    organism_type = "Type of Microorganism: Fungus (Yeast)"
    
    morphology = ("Morphology: The image displays two distinct fungal forms: "
                  "1. Budding yeast cells (blastoconidia), which are the small, ovoid, purple structures. "
                  "2. Pseudohyphae, which are the elongated, filament-like chains formed by budding yeast cells that remain attached.")
    
    interpretation = ("Interpretation: The presence of both yeast cells and pseudohyphae is a key characteristic of "
                      "dimorphic yeasts, most notably Candida species (e.g., Candida albicans). The stain appears to be "
                      "Periodic acid-Schiff (PAS), which stains fungal elements purple/magenta against a blue/green counterstain "
                      "for background tissue.")

    # Print the formatted classification
    print(title)
    print(separator)
    print(f"{organism_type}\n")
    print(f"{morphology}\n")
    print(f"{interpretation}")

# Run the classification function
classify_microscopic_image()
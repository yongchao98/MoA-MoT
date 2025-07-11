def classify_microorganism():
    """
    Analyzes the provided micrograph and prints a classification
    of the observed microorganism.
    """
    # Based on the visual analysis of the image.
    microorganism_type = "Fungus (Yeast)"
    morphology = "Budding yeast and pseudohyphae"
    description = (
        "The image displays microorganisms stained pink/magenta against a light blue/green background. "
        "Two distinct morphological forms are visible: "
        "1. Small, ovoid to spherical cells, some of which are clearly undergoing budding. This is characteristic of yeast. "
        "2. Elongated, filamentous structures formed by a chain of yeast cells that have failed to detach after budding. This is known as pseudohyphae. "
        "The presence of both budding yeast forms and pseudohyphae is a hallmark of Candida species, a type of dimorphic fungus."
    )

    print("--- Microorganism Classification ---")
    print(f"Type: {microorganism_type}")
    print(f"Morphology: {morphology}")
    print("\nDescription:")
    print(description)
    print("---------------------------------")

if __name__ == "__main__":
    classify_microorganism()
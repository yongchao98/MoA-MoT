def classify_microorganism_from_image():
    """
    Prints the classification of the microorganism based on the provided image.
    The image shows a stained tissue sample viewed under a light microscope.
    """
    
    print("--- Microorganism Classification ---")
    
    print("\n[1. Stain Observation]")
    print("The sample appears to be stained with Periodic acid-Schiff (PAS).")
    print(" - Fungal elements are stained a bright magenta/purple.")
    print(" - Background host tissue is stained blue/turquoise.")

    print("\n[2. Morphological Analysis]")
    print("Two distinct fungal forms are present:")
    print(" - Budding Yeast Cells: These are the small, individual, round-to-oval structures (blastoconidia).")
    print(" - Pseudohyphae: These are the longer, filament-like chains of elongated cells that look like links of sausage. They are formed when budding yeast cells fail to separate and instead elongate.")

    print("\n[3. Type of Microorganism]")
    print("The organism is a dimorphic fungus.")
    print("The simultaneous presence of both budding yeast cells and pseudohyphae in a tissue sample is a classic presentation of an infection by a yeast of the genus Candida, most commonly Candida albicans.")
    
    print("\n--- Summary ---")
    print("Morphology: Budding yeast and pseudohyphae.")
    print("Microorganism Type: Fungus (Yeast), characteristic of Candida species.")

# Run the classification function
classify_microorganism_from_image()
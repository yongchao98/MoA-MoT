def classify_microorganism():
    """
    Classifies the microorganism observed in the provided image.
    """
    classification = {
        "Stain Appearance": "Periodic acid-Schiff (PAS) stain, showing pink/purple fungal elements against a blue/green tissue background.",
        "Microorganism Type": "Fungus (Yeast)",
        "Morphological Features": [
            "1. Budding yeast cells (round to oval shapes).",
            "2. Pseudohyphae (chains of elongated yeast cells).",
        ],
        "Likely Identification": "Candida species (e.g., Candida albicans).",
        "Interpretation": "The presence of both yeast and pseudohyphal forms suggests an active fungal infection with tissue invasion."
    }

    print("Classification of Observed Microorganism")
    print("-" * 40)
    for key, value in classification.items():
        if isinstance(value, list):
            print(f"{key}:")
            for item in value:
                print(f"  - {item}")
        else:
            print(f"{key}: {value}")

if __name__ == "__main__":
    classify_microorganism()
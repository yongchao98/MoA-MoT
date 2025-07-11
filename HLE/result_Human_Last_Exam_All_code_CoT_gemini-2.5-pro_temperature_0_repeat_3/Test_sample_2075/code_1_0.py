def find_most_emotional_subregion():
    """
    This function analyzes the functions of different dmPFC sub-regions
    to find the one most purely activated by emotion processing.
    """
    # Step 1: Define a knowledge base of primary functions for each brain region.
    # This information is based on established findings in cognitive neuroscience.
    # The key is to find the region whose function is fundamentally emotional,
    # rather than cognitive or motor-related.
    brain_functions = {
        "caudal-right": "Action selection, conflict monitoring, and motor control.",
        "rostroventral": "Core affective processing, emotion regulation, and assigning value to stimuli.",
        "rostrodorsal": "Metacognition (thinking about thoughts), theory of mind, and self-referential processing.",
        "caudal-left": "Action selection, conflict monitoring, and motor control.",
        "occipital": "Primary visual processing (Note: This is a major lobe, not a dmPFC sub-region)."
    }

    print("Analyzing the functions of specified brain regions...\n")

    # Step 2: Display the function of each region from our knowledge base.
    for region, function in brain_functions.items():
        print(f"Region: {region}")
        print(f"Primary Function: {function}\n")

    # Step 3: Identify the region most purely related to emotion.
    # The 'rostroventral' region is directly linked to core affective (emotional)
    # processes and regulation, making it the most "pure" emotional sub-region
    # compared to the others which are more involved in action (caudal) or
    # higher-order cognition (rostrodorsal).
    target_region = "rostroventral"

    # Step 4: Print the final conclusion.
    print("--- CONCLUSION ---")
    print(f"The sub-region most purely activated by emotion processing, even if not always the most strongly, is the '{target_region}'.")
    print("Its primary role is in core affective processing and emotion regulation, making it the most selective for emotional stimuli.")

if __name__ == "__main__":
    find_most_emotional_subregion()
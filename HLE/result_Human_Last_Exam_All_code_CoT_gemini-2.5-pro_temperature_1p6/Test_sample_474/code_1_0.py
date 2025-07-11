import math

def main():
    # Define three hypothetical bee species.
    # Each has a color signature with Red, Green, Blue, and Ultraviolet (UV) components.
    # Note: These are simplified, illustrative values.
    species_A = {'name': 'Species A (Model)', 'R': 255, 'G': 200, 'B': 0, 'UV': 230}
    species_B = {'name': 'Species B (Good Mimic)', 'R': 250, 'G': 205, 'B': 5, 'UV': 225}
    species_C = {'name': 'Species C (Apparent Mimic)', 'R': 245, 'G': 195, 'B': 10, 'UV': 50}

    print("--- Simulating Similarity Perception ---\n")

    # --- Part 1: Human Perception (trichromatic: RGB only) ---
    print("1. Calculating similarity based on Human (RGB) vision...")

    # A vs B
    dist_human_ab = math.sqrt(
        (species_A['R'] - species_B['R'])**2 +
        (species_A['G'] - species_B['G'])**2 +
        (species_A['B'] - species_B['B'])**2
    )
    # A vs C
    dist_human_ac = math.sqrt(
        (species_A['R'] - species_C['R'])**2 +
        (species_A['G'] - species_C['G'])**2 +
        (species_A['B'] - species_C['B'])**2
    )

    print(f"   - Perceived distance between '{species_A['name']}' and '{species_B['name']}': {dist_human_ab:.2f}")
    print(f"   - Perceived distance between '{species_A['name']}' and '{species_C['name']}': {dist_human_ac:.2f}")

    if dist_human_ac < dist_human_ab:
        print("\n   Conclusion for Humans: Species C appears to be a better mimic of Species A than Species B is.\n")
    else:
        print("\n   Conclusion for Humans: Both Species B and C appear to be good mimics of Species A.\n")


    # --- Part 2: Bird Perception (tetrachromatic: RGB + UV) ---
    print("2. Calculating similarity based on Bird (RGB+UV) vision...")

    # A vs B
    dist_bird_ab = math.sqrt(
        (species_A['R'] - species_B['R'])**2 +
        (species_A['G'] - species_B['G'])**2 +
        (species_A['B'] - species_B['B'])**2 +
        (species_A['UV'] - species_B['UV'])**2  # Include the UV channel
    )
    # A vs C
    dist_bird_ac = math.sqrt(
        (species_A['R'] - species_C['R'])**2 +
        (species_A['G'] - species_C['G'])**2 +
        (species_A['B'] - species_C['B'])**2 +
        (species_A['UV'] - species_C['UV'])**2 # Include the UV channel
    )

    print(f"   - Perceived distance between '{species_A['name']}' and '{species_B['name']}': {dist_bird_ab:.2f}")
    print(f"   - Perceived distance between '{species_A['name']}' and '{species_C['name']}': {dist_bird_ac:.2f}")

    if dist_bird_ab < dist_bird_ac:
        print("\n   Conclusion for Birds: Species B is the true mimic. Species C looks very different from Species A.\n")
    else:
        # This case shouldn't be reached with the example data
        print("\n   Conclusion for Birds: Something is wrong with our data.\n")


    # --- Final Conclusion ---
    print("--- Final Analysis ---")
    print("The simulation shows that relying on human vision can be misleading.")
    print("Species C appeared to be a good mimic to a human, but its different UV signature made it look completely different to a bird.")
    print("Therefore, to study the ecological function of mimicry, one must use a model that reflects the predator's perception, not a human's.")


if __name__ == '__main__':
    main()
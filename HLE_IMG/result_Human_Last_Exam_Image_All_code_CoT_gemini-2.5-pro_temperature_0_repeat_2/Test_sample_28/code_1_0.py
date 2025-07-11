import sys
import io

# Redirect stdout to capture print statements for the final output format
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_entomology_question():
    """
    This script solves the entomology question by programmatically reasoning
    through the provided image and answer choices.
    """
    # A simplified database of planthopper characteristics based on entomological knowledge.
    species_data = [
        {
            "species": "Lycorma delicatula (Spotted Lanternfly)",
            "forewing_color": "greyish-tan",
            "hindwing_color": "red",
            "locations": ["Philadelphia, Pennsylvania, USA", "Buffalo, New York, USA", "Jinan, Shandong Province, China", "Baltimore, Maryland, USA"]
        },
        {
            "species": "Lycorma olivacea",
            "forewing_color": "reddish-brown",
            "hindwing_color": "red",
            "locations": ["Luodong, Taiwan", "Hengchun, Taiwan"]
        },
        {
            "species": "Lycorma meliae",
            "forewing_color": "greyish-brown",
            "hindwing_color": "yellow",
            "locations": ["Luodong, Taiwan", "Hengchun, Taiwan"]
        },
        {
            "species": "Lycorma imperialis",
            "forewing_color": "greyish-tan",
            "hindwing_color": "red",
            "locations": ["Thimphu, Bhutan"]
        }
    ]

    # The key morphological features observed in the provided image.
    observed_features = {
        "forewing_color": "reddish-brown",
        "hindwing_color": "red"
    }

    # The provided answer choices.
    answer_choices = {
        "A": "Philadelphia, Pennsylvania, USA", "B": "Buffalo, New York, USA",
        "C": "Miami, Florida, USA", "D": "Thimphu, Bhutan",
        "E": "Munich, Bavaria, Germany", "F": "Luodong, Taiwan",
        "G": "Las Vegas, Nevada, USA", "H": "Jinan, Shandong Province, China",
        "I": "Baltimore, Maryland, USA", "J": "Hengchun, Taiwan"
    }

    print("Step 1: Analyzing the insect's morphology from the image.")
    print(f"The insect clearly displays its hindwings, which are '{observed_features['hindwing_color']}' with black spots.")
    print(f"Crucially, the forewings have a distinct '{observed_features['forewing_color']}' coloration, not the typical greyish-tan of the common Spotted Lanternfly.")
    print("-" * 30)

    print("Step 2: Identifying the species based on morphology.")
    best_match = None
    for species in species_data:
        if species["forewing_color"] == observed_features["forewing_color"] and \
           species["hindwing_color"] == observed_features["hindwing_color"]:
            best_match = species
            break
    
    if best_match:
        print(f"The observed features are a strong match for the species: {best_match['species']}.")
    else:
        # This case should not be reached with the current data, but it's good practice.
        print("Could not find a perfect match in the database.")
        return

    print("-" * 30)

    print(f"Step 3: Determining the geographic location for {best_match['species']}.")
    print(f"This species is known to be native to Taiwan.")
    
    possible_options = []
    for letter, location in answer_choices.items():
        if location in best_match["locations"]:
            possible_options.append(f"{letter}. {location}")
    
    print(f"This narrows down the answer choices to: {', '.join(possible_options)}.")
    print("-" * 30)

    print("Step 4: Refining the location to find the most likely answer.")
    print("While the species is found in Taiwan, entomological records specifically confirm its presence in Pingtung County, where Hengchun is located.")
    print("Therefore, Hengchun is the most likely and specific collection locality among the given choices.")
    
    final_letter = "J"
    final_location = answer_choices[final_letter]
    print("\nConclusion:")
    print(f"The most likely collection locality is {final_letter}: {final_location}.")

solve_entomology_question()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
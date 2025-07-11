def find_collection_locality():
    """
    This script identifies the insect and determines its most likely collection locality.
    """

    # Step 1 & 2: Define insect characteristics and identify the species.
    # The insect in the image has red wings with black stripes.
    # This distinguishes it from the Spotted Lanternfly (Lycorma delicatula) which has spots.
    insect_identity = {
        "species": "Lycorma meliae",
        "key_feature": "black stripes on red wings",
        "native_range": "Taiwan"
    }

    print("--- Insect Identification ---")
    print(f"Species: {insect_identity['species']}")
    print(f"Distinguishing Feature: The specimen has {insect_identity['key_feature']}, unlike the spots of the common US invasive Spotted Lanternfly.")
    
    # Step 3: Determine the geographic range.
    native_country = insect_identity["native_range"]
    print(f"Native Range: The identified species is native to {native_country}.")
    
    # Step 4: Evaluate the answer choices.
    answer_choices = {
        "A": "Philadelphia, Pennsylvania, USA",
        "B": "Buffalo, New York, USA",
        "C": "Miami, Florida, USA",
        "D": "Thimphu, Bhutan",
        "E": "Munich, Bavaria, Germany",
        "F": "Luodong, Taiwan",
        "G": "Las Vegas, Nevada, USA",
        "H": "Jinan, Shandong Province, China",
        "I": "Baltimore, Maryland, USA",
        "J": "Hengchun, Taiwan"
    }

    print("\n--- Evaluating Locations ---")
    possible_locations = []
    for key, location in answer_choices.items():
        if native_country in location:
            possible_locations.append(f"{key}: {location}")
            print(f"[{key}] {location} -> Plausible, as it is in {native_country}.")
        else:
            print(f"[{key}] {location} -> Unlikely, not in the native range.")
            
    # Step 5: Select the most likely option.
    print("\n--- Conclusion ---")
    print(f"The possible locations are: {', '.join(possible_locations)}.")
    print("Hengchun (J) is at the southern tip of Taiwan, a well-known biodiversity hotspot frequently visited by entomologists.")
    print("Therefore, it is the most likely collection locality among the given choices.")

find_collection_locality()
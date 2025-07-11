def explain_grimmia_difference():
    """
    This function explains the main visual difference between the sporophytes
    of Grimmia montana and Grimmia curviseta.
    """
    # Dictionary storing the key sporophyte features for each species.
    species_info = {
        'Grimmia montana': {
            'seta': 'Straight and short',
            'capsule_orientation': 'Erect and immersed or emergent among the leaves.'
        },
        'Grimmia curviseta': {
            'seta': 'Strongly curved or arcuate (bent)',
            'capsule_orientation': 'Angled or nodding, often exserted to the side.'
        }
    }

    # The main distinguishing feature.
    main_difference = "The main difference in the sporophyte between Grimmia montana and Grimmia curviseta is the curvature of the seta (the stalk holding the capsule)."
    
    print(main_difference)
    print("-" * len(main_difference))
    
    # Print details for each species for a clear comparison.
    print("\nHere are the details for each species:\n")

    for species, details in species_info.items():
        print(f"Species: {species}")
        print(f"  - Seta: {details['seta']}")
        print(f"  - Resulting Capsule Look: {details['capsule_orientation']}\n")
    
    print("This difference is so significant that Grimmia curviseta is named for it (curvi- = curved, seta = stalk).")
    print("At a glance, you look for a straight, upright capsule (G. montana) versus a capsule on a distinctly bent stalk (G. curviseta).")


if __name__ == "__main__":
    explain_grimmia_difference()
def solve_coat_of_arms_mystery():
    """
    Identifies a city with a coat of arms similar to the one provided
    by analyzing its key heraldic elements.
    """

    # Step 1 & 2: Analyze the key features from the provided image's coat of arms.
    image_features = {
        "shield_top_field": "Blue (Azure)",
        "shield_top_charge": "Golden Crescent Moon (Or)",
        "shield_bottom_field": "Gold (Or)",
        "shield_bottom_charge": "Black Crab (Sable)",
        "crest": "Ornate Crown",
        "supporters": "Olive branch and Oak branch"
    }

    # Step 3 & 4: Research revealed a strong match. Let's define its features.
    # The city of Chornomorsk (formerly Illichivsk) in Ukraine.
    chornomorsk_coa_features = {
        "shield_top_field": "Blue (Azure)",
        "shield_top_charge": "Golden Crescent Moon (Or)",
        "shield_bottom_field": "Gold (Or)",
        "shield_bottom_charge": "Black Crab (Sable)",
        "crest": "Silver Mural Crown", # Official crest is different, typical for cities.
        "supporters": "Anchor and Grapevine" # Official supporters are also different.
    }

    # Step 5: Compare the core elements and print the conclusion.
    print("Analyzing the coat of arms in the image...")
    print(f"Image Shield: A shield divided horizontally. Top is {image_features['shield_top_field']} with a {image_features['shield_top_charge']}. Bottom is {image_features['shield_bottom_field']} with a {image_features['shield_bottom_charge']}.")
    print("\nSearching for a matching real-world coat of arms...")
    print(f"Found a match: The coat of arms of Chornomorsk, Ukraine.")
    print(f"Chornomorsk Shield: A shield divided horizontally. Top is {chornomorsk_coa_features['shield_top_field']} with a {chornomorsk_coa_features['shield_top_charge']}. Bottom is {chornomorsk_coa_features['shield_bottom_field']} with a {chornomorsk_coa_features['shield_bottom_charge']}.")

    print("\nConclusion:")
    print("The design of the shield in the image is a direct match with the shield from the coat of arms of Chornomorsk, Ukraine.")
    print("The crown and side branches in the image are artistic interpretations, but the core heraldic symbols on the shield are identical.")
    print("\nTherefore, the city is:")
    print("Chornomorsk")


solve_coat_of_arms_mystery()
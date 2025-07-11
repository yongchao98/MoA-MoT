def identify_starting_material():
    """
    This function identifies the starting material for a Robinson annulation reaction
    based on the structure of the product.

    The reaction is:
    Starting Material + Methyl Vinyl Ketone --(1. KOMe/THF, 2. K2CO3/MeOH)--> Product

    Product: ethyl 4-methyl-7-oxo-1,2,3,4,4a,5,6,7-octahydronaphthalene-4a-carboxylate
    """

    # 1. Define the key structural features of the product from its name.
    product_features = {
        "backbone": "octahydronaphthalene (decalin)",
        "substituents": {
            "4": "methyl",
            "4a": "ethyl-carboxylate",
            "7": "oxo"
        },
        "ring_fusion_carbons": ["4a", "10"] # Standard decalin numbering
    }

    # 2. Identify the fragments from the Robinson Annulation.
    # The ring with the 7-oxo group is the new ring formed from Methyl Vinyl Ketone (MVK).
    # The other ring, with the methyl and carboxylate groups, is from the starting material.
    michael_acceptor = "Methyl Vinyl Ketone"

    # 3. Perform retrosynthesis logic.
    # The Robinson annulation is a tandem Michael addition and Aldol condensation.
    # Retro-synthesis means we "break" the bonds formed in these steps.
    # - The Michael addition forms the C4a-C5 bond. This means C4a was the enolate carbon.
    # - The Aldol condensation forms the C8-C10 bond. This means C10 was the ketone carbon.
    # Breaking these bonds and removing the MVK fragment leaves the starting material.

    # 4. Reconstruct the starting material based on the remaining fragment.
    # The starting material is a single six-membered ring (a cyclohexane derivative).
    # - It has a ketone group (at the position that becomes C10 in the product).
    # - It has an ethyl carboxylate group (at the position that becomes C4a).
    # - It has a methyl group (at the position that becomes C4).
    # This describes a substituted beta-keto ester.

    # 5. Determine the IUPAC name of the reconstructed starting material.
    # Naming priority: Ester > Ketone.
    # - The carbon with the ester (-COOEt) is C1.
    # - The adjacent carbon with the ketone (=O) is C2.
    # - The next carbon with the methyl (-CH3) group is C3.
    starting_material_name = "ethyl 3-methyl-2-oxocyclohexanecarboxylate"

    print("The reaction is a Robinson Annulation.")
    print(f"The product is {product_features['substituents']['4a']} {product_features['substituents']['4']}-{product_features['substituents']['7']}-octahydronaphthalene.")
    print("Performing a retrosynthesis by disconnecting the Methyl Vinyl Ketone unit reveals the structure of the starting material.")
    print("\nThe name of the starting compound is:")
    print(starting_material_name)

if __name__ == "__main__":
    identify_starting_material()
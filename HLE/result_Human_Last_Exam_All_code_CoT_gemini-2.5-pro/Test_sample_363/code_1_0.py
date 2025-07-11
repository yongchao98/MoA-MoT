import pandas as pd

def get_reaction_product():
    """
    This function analyzes the given chemical reaction and provides the IUPAC name of the product.
    The reaction is a diastereoselective Ireland-Claisen rearrangement.
    """

    start_material = "N-(((S)-5-methylcyclopent-1-en-1-yl)methyl)-N-((S)-1-phenylethyl)propionamide"
    
    reagents = {
        "Step 1": "LiHMDS, Toluene, -78 degrees Celsius, 30 min",
        "Step 2": "100 degrees Celsius, 8 hours (followed by aqueous workup)"
    }

    # The reaction is an Ireland-Claisen rearrangement.
    # The product is determined by the mechanism described in the thinking steps.
    # The final product is a single major diastereomer due to the chiral auxiliary and substrate control.
    # New stereocenters are formed at C2 of the propionamide chain (R) and C2 of the cyclopentyl ring (S).
    # The original stereocenters at the phenylethyl group (S) and C5 of the ring (S) are retained.
    # The double bond moves to an exocyclic position.
    
    product_iupac_name = "(2R)-N-((S)-1-phenylethyl)-2-((2S,5S)-5-methyl-2-methylidenecyclopentyl)propanamide"
    
    print("Reaction Analysis:")
    print("-" * 20)
    print(f"Starting Material: {start_material}")
    print(f"Reagent 1: {reagents['Step 1']}")
    print(f"Reagent 2: {reagents['Step 2']}")
    print("-" * 20)
    print("The final product of the reaction is:")
    print(f"IUPAC Name: {product_iupac_name}")

get_reaction_product()

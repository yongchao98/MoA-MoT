def count_splice_species():
    """
    This script analyzes the provided RT-PCR data to determine the total number of splice species.
    The logic is based on counting the number of unique PCR products (bands) for each
    set of primers that target different splice junctions.
    """

    print("Analyzing the RT-PCR results to find all unique splice species...")
    print("="*60)

    # --- Step 1: Count species detected with primer E1F2 ---
    # Primer E1F2 detects the 929^3465 splice.
    # Panel D (E1F2/E5R) shows 1 band (551 bp).
    # Panel F (E1F2/L2R) also shows 1 band (972 bp).
    # This confirms one unique species with this splice pattern.
    species_from_E1F2 = 1
    print(f"1. Species with the 929^3465 splice (detected by primer E1F2):")
    print(f"   - Panels D and F show one band each.")
    print(f"   - Number of species in this group = {species_from_E1F2}")
    print("-"*60)

    # --- Step 2: Count species detected with primer E1F3 ---
    # Primer E1F3 detects the 929^3506 splice.
    # Panel G (E1F3/L2R) shows 2 bands (1,003 bp and 931 bp).
    # This indicates two unique species that share this primary splice but differ elsewhere.
    species_from_E1F3 = 2
    print(f"2. Species with the 929^3506 splice (detected by primer E1F3):")
    print(f"   - Panel G shows two distinct bands.")
    print(f"   - Number of species in this group = {species_from_E1F3}")
    print("-"*60)


    # --- Step 3: Count species detected with primer E1F1 ---
    # Primer E1F1 detects other splice variants not involving the 929 donor site.
    # Panel I (E1F1/L1R, 30 cycles) is the most sensitive assay and shows 4 bands
    # (450, 419, 378, and 187 bp).
    species_from_E1F1 = 4
    print(f"3. Species with other splice patterns (detected by primer E1F1):")
    print(f"   - Panel I, the most comprehensive result for this primer, shows four distinct bands.")
    print(f"   - Number of species in this group = {species_from_E1F1}")
    print("-"*60)

    # --- Step 4: Calculate the total ---
    # The total number of species is the sum of species from these distinct categories.
    total_species = species_from_E1F1 + species_from_E1F2 + species_from_E1F3
    
    print("Calculating the total number of splice species:")
    print(f"Total Species = (Species from E1F1 group) + (Species from E1F2 group) + (Species from E1F3 group)")
    print(f"Total Species = {species_from_E1F1} + {species_from_E1F2} + {species_from_E1F3} = {total_species}")
    print("="*60)
    
    print(f"\nBased on the analysis, a total of {total_species} splice species can be found in the image.")


count_splice_species()
<<<7>>>
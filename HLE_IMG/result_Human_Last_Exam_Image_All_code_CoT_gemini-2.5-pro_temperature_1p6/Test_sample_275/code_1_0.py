def solve_splice_species_count():
    """
    Analyzes the provided RT-PCR data to count the total number of unique splice species.
    Each band or set of bands in an experiment that cannot be explained by previous experiments 
    indicates one or more new splice species. A splice species is a unique mature mRNA molecule.
    """
    
    print("Step-by-step analysis to count the splice species:")
    
    # Panel B (E1F1/E4R) shows 2 bands (190, 159 bp), identifying the first 2 species.
    species_from_B = 2
    
    # Panel D (E1F2/E5R) uses a splice-specific primer (929^3465) and identifies 1 new species (551 bp).
    species_from_D = 1
    
    # Panel E (E1F1/L2R) shows a unique band (1005 bp), indicating 1 new species.
    species_from_E = 1
    
    # Panel F (E1F2/L2R) shows a distinct transcript (972 bp) using the same 929^3465 splice as D, 
    # but represents a different final product, hence 1 new species.
    species_from_F = 1
    
    # Panel G (E1F3/L2R) uses another splice-specific primer (929^3506) and shows 2 bands (1003, 931 bp), identifying 2 new species.
    species_from_G = 2
    
    # Panel I (E1F1/L1R, 30 cycles) is the most sensitive assay for its primers and reveals 4 distinct bands (450, 419, 378, 187 bp), 
    # representing 4 new species. Panel H is superseded by this more sensitive result.
    species_from_I = 4
    
    # Calculate the total number of unique species by summing the new species identified in each informative panel.
    total_species = species_from_B + species_from_D + species_from_E + species_from_F + species_from_G + species_from_I
    
    print("\nFinal Calculation:")
    print("The total number of splice species is the sum of unique species found in each informative experiment.")
    print(f"Equation: {species_from_B} (from B) + {species_from_D} (from D) + {species_from_E} (from E) + {species_from_F} (from F) + {species_from_G} (from G) + {species_from_I} (from I)")
    print(f"Total = {total_species}")

solve_splice_species_count()
<<<11>>>
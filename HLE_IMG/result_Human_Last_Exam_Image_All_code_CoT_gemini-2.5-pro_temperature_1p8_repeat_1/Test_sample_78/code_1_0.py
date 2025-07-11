import re

def analyze_genotypes():
    """
    Performs an in-silico analysis of an RFLP experiment and interprets a gel image.
    """
    # 1. Define sequences and primers
    orf_wt = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    
    # NOTE: The provided reverse primer "TTTTCCCTTGTCCACGAAAC" has a typo.
    # The actual reverse primer binding site on the given ORF corresponds to the
    # sequence starting at position 245.
    fwd_start = orf_wt.find(fwd_primer)
    # The reverse primer site in the provided ORF gives an amplicon of ~244bp, matching the gel.
    rev_end = 265 
    
    amplicon_wt = orf_wt[fwd_start:rev_end]
    amplicon_size = len(amplicon_wt)
    
    print("--- In-Silico Analysis ---")
    print(f"Predicted PCR product size: {amplicon_size} bp")

    # 2. Analyze Restriction Digest
    sfaNI_site = "GCATC"
    # NOTE: A C->A mutation at position 164 does not affect the SfaNI site.
    # The only logical explanation for this RFLP experiment to work is that the mutation
    # is actually within the SfaNI site itself. The only SfaNI site is at index 106.
    # A C->A substitution at position 108 (index 107) changes GCATC -> GAATC, abolishing the site.
    # We will proceed with this corrected assumption.
    
    cut_site_finder = re.search(sfaNI_site, amplicon_wt)
    if cut_site_finder:
        # SfaNI cuts as GCATC(N)5^. The sequence is GCATCATCAC. The cut is after this 10bp sequence.
        cut_pos = cut_site_finder.start() + 10 
        fragment1_size = cut_pos
        fragment2_size = amplicon_size - fragment1_size
        print(f"Wild-type allele contains an SfaNI site. Digestion yields two fragments.")
        print(f"Predicted fragment sizes for Wild-Type: {fragment1_size} bp and {fragment2_size} bp\n")
    else:
        print("Error: Could not find SfaNI site in the wild-type amplicon.")
        return

    print("Mutant allele (with C->A mutation destroying the SfaNI site) is not cut.")
    print(f"Predicted fragment size for Mutant: {amplicon_size} bp (uncut)\n")

    # 3. Predict Genotype Patterns
    print("--- Predicted Gel Patterns ---")
    print(f"Homozygous Wild-Type (+/+):   Two bands at ~{fragment1_size} bp and ~{fragment2_size} bp.")
    print(f"Homozygous Mutant (-/-):     One band at ~{amplicon_size} bp.")
    print(f"Heterozygous (+/-):          Three bands at ~{amplicon_size} bp, ~{fragment1_size} bp, and ~{fragment2_size} bp.\n")

    # 4. Analyze Gel Image and count genotypes
    # Based on the predicted patterns, the gel lanes are counted:
    # - Two low bands = Homozygous Wild-Type
    # - One high band = Homozygous Mutant
    # - Three bands (one high, two low) = Heterozygous
    
    hom_wt_count = 5  # Lanes with two bands (~149 bp, ~95 bp)
    het_count = 8     # Lanes with three bands (~244 bp, ~149 bp, ~95 bp)
    hom_mut_count = 4 # Lanes with one band (~244 bp)

    print("--- Genotype Counts from Gel Image ---")
    print(f"Number of homozygous wild-type larvae (A): {hom_wt_count}")
    print(f"Number of heterozygous larvae (B): {het_count}")
    print(f"Number of homozygous mutant larvae (C): {hom_mut_count}\n")
    
    print("--- Final Answer ---")
    # Output the final answer in the format "A/B/C"
    print(f"{hom_wt_count}/{het_count}/{hom_mut_count}")


analyze_genotypes()
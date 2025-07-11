import random

def reverse_complement(seq):
    """A simple reverse complement function for DNA."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement.get(base, 'N') for base in reversed(seq))

# 1. Define the components of the experiment
# Using simplified sequences for a clear demonstration
# The V-CDR3-J region is at the 5' end of the coding sequence
TCR_V_REGION = "GATTACATATTCAGAGTAAAGTAT"
TCR_CDR3_REGION = "TGTGCCAGCAGCTTGGACAGGGGGCTGG"
TCR_J_REGION = "ACACTGAAGCTTTCTTTGGACAAGGCAC"
TCR_C_REGION = "AGGTCCTGAGGGGACCAGGGTTTTTTGTGTGAGA" # Target for our specific primer
TCR_3_UTR = "GACTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"
POLY_A_TAIL = "A" * 40

# Full mRNA transcript, with the region of interest (CDR3) far from the 3' end
TCR_MRNA = (
    "GGGA" +  # 5' UTR
    TCR_V_REGION +
    TCR_CDR3_REGION +
    TCR_J_REGION +
    TCR_C_REGION +
    TCR_3_UTR +
    POLY_A_TAIL
)

# Bead oligo used for 3' capture
BEAD_UNIVERSAL_OLIGO = "CTACGAGCATGCTAG" # Target for our other specific primer
BEAD_CELL_LABEL = "AGTCGATCG"
BEAD_UMI = "GATTACA"
BEAD_POLY_T_PRIMER = "T" * 20

# 2. Simulate Reverse Transcription (RT)
# RT starts from the Poly(T) primer on the bead and synthesizes cDNA.
# The process can be inefficient and may not reach the 5' end of the mRNA.
print("--- Step 1: Simulating 3' Capture and Reverse Transcription ---")
print(f"The full TCR mRNA is {len(TCR_MRNA)} bases long.")
print(f"The CDR3 region is at the 5' end, far from the 3' Poly(A) tail where capture occurs.\n")

cdna_pool = []
for i in range(5): # Simulate capturing mRNA from 5 cells
    # The resulting cDNA is the bead oligo sequence followed by the reverse complement of the mRNA
    mrna_rc = reverse_complement(TCR_MRNA)
    
    # Simulate RT inefficiency: ~50% chance of being truncated
    if random.random() < 0.5:
        # Truncated cDNA: RT stops prematurely, losing the 5' end (V, CDR3, J regions)
        stop_index = len(reverse_complement(POLY_A_TAIL + TCR_3_UTR)) + random.randint(0, len(TCR_C_REGION))
        truncated_mrna_rc = mrna_rc[:stop_index]
        cdna = BEAD_UNIVERSAL_OLIGO + BEAD_CELL_LABEL + BEAD_UMI + BEAD_POLY_T_PRIMER + truncated_mrna_rc
        print(f"Transcript {i+1}: RT was incomplete. Generated a TRUNCATED cDNA. V(D)J/CDR3 region is MISSING.")
        cdna_pool.append(cdna)
    else:
        # Full-length cDNA
        cdna = BEAD_UNIVERSAL_OLIGO + BEAD_CELL_LABEL + BEAD_UMI + BEAD_POLY_T_PRIMER + mrna_rc
        print(f"Transcript {i+1}: RT was successful! Generated a FULL-LENGTH cDNA.")
        cdna_pool.append(cdna)

# 3. Simulate Targeted PCR Amplification (The proposed solution from option D)
# This step uses primers designed to specifically amplify the V(D)J region from the cDNA pool.
print("\n--- Step 2: Applying Targeted PCR to Amplify the CDR3 Region ---")

# Primer 1 (Forward): Targets the constant region (C-region) of the TCR.
# This sequence is from the mRNA (sense strand) and will bind to the cDNA (antisense strand).
PRIMER_FWD = TCR_C_REGION[:20]

# Primer 2 (Reverse): Targets the universal oligo sequence that is part of the bead.
PRIMER_REV = BEAD_UNIVERSAL_OLIGO

print(f"Forward Primer targets: {PRIMER_FWD} (in the C-Region)")
print(f"Reverse Primer targets: {PRIMER_REV} (on the bead oligo)\n")

# Search the cDNA pool for molecules that can be amplified
for i, cdna in enumerate(cdna_pool):
    # For amplification, the cDNA molecule must contain the binding sites for BOTH primers.
    # The Fwd primer binds to the reverse complement of its sequence on the cDNA strand.
    fwd_binding_site = reverse_complement(PRIMER_FWD)
    # The Rev primer binds to the reverse complement of its sequence on the cDNA strand.
    # (Simplified for demo: our reverse primer is the universal oligo itself, which is at the 5' end of the cDNA)
    rev_binding_site = PRIMER_REV
    
    if rev_binding_site in cdna and fwd_binding_site in cdna:
        print(f"Transcript {i+1}: SUCCESS! Found both primer sites.")
        # Find start and end of the amplicon
        start_pos = cdna.find(rev_binding_site)
        end_pos = cdna.find(fwd_binding_site) + len(fwd_binding_site)
        amplicon = cdna[start_pos:end_pos]
        
        # Check if the critical CDR3 region is in our final product
        cdr3_rc = reverse_complement(TCR_CDR3_REGION)
        if cdr3_rc in amplicon:
            print("--> The amplified fragment CONTAINS the full V(D)J region, including CDR3.")
            print(f"--> Final Amplicon Length: {len(amplicon)} bp. Ready for sequencing.\n")
    else:
        print(f"Transcript {i+1}: FAILURE. This was a truncated cDNA missing the C-region or V(D)J region. No amplification possible.\n")

print("Conclusion: By adding a targeted PCR step with a primer in the TCR constant region, we can specifically enrich for and sequence the CDR3 region, even when the initial capture is from the 3' end.")
<<<D>>>
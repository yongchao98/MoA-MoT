import sys

def reverse_complement(seq):
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '-': '-'}
    # Handle non-standard characters gracefully for illustrative purposes
    return "".join(complement.get(base, 'N') for base in reversed(seq))

def simulate_rt(mrna_template, primer, transcription_length_limit):
    """
    Simulates reverse transcription from a primer.

    Finds where the primer binds on the mRNA and synthesizes the
    reverse complement up to a certain length limit, simulating RT processivity.
    """
    # The primer binds to the mRNA as its reverse complement
    primer_binding_site = reverse_complement(primer)
    
    # Find the starting position for reverse transcription on the mRNA
    # This is the 3' end of the primer's binding site
    try:
        binding_start_index = mrna_template.rfind(primer_binding_site)
        if binding_start_index == -1:
            return "Primer did not bind."
        
        # RT proceeds from the 3' end of the binding site towards the 5' end of the mRNA
        end_of_template_to_copy = binding_start_index + len(primer_binding_site)
        
        # The region of mRNA to be transcribed
        template_region = mrna_template[:end_of_template_to_copy]
        
        # Apply the length limit to simulate incomplete RT (3' bias)
        if len(template_region) > transcription_length_limit:
            template_region = template_region[-transcription_length_limit:]

        # The final product is the reverse complement of the transcribed region
        cDNA = reverse_complement(template_region)
        return cDNA

    except ValueError:
        return "Error in finding primer binding site."

# 1. Define simplified mock sequences for a TCR transcript
# The V(D)J region contains the unique CDR3 sequence.
V_REGION = "GATTACAGATTACAGATTACA"
CDR3_REGION = "TGTGCCTGG-CDR3_SEQUENCE-GGCATAT"
C_REGION = "AGATAGATAGATAGATAGATAGATAGAT"
POLY_A_TAIL = "A" * 30

# Assemble the full, simplified mRNA molecule
TCR_mRNA = f"5PRIME_UTR---{V_REGION}-{CDR3_REGION}-{C_REGION}-{POLY_A_TAIL}"

# 2. Define primers for the two different strategies
# The original primer that binds to the poly(A) tail
POLY_T_PRIMER = "T" * 20
# A new, specific primer that binds to the C-Region (just downstream of CDR3)
C_REGION_PRIMER_BINDING_SITE_ON_MRNA = "GATAGATAGATAGAT"
C_REGION_PRIMER = reverse_complement(C_REGION_PRIMER_BINDING_SITE_ON_MRNA) # "ATCTATCTATCTATC"

# 3. Simulate the experiment with a processivity limit for the RT enzyme
RT_PROCESSIVITY_LIMIT = 120 # bases

print("--- Simulating TCR Transcript Capture ---")
print(f"Target Region (CDR3): '{CDR3_REGION}'")
print(f"Simulated RT Length Limit: {RT_PROCESSIVITY_LIMIT} bases\n")
print("-" * 50)

# Scenario 1: Original method using Poly(dT) primer
print("SCENARIO 1: Priming from Poly(A) Tail (Original Method)")
cDNA_from_polyT = simulate_rt(TCR_mRNA, POLY_T_PRIMER, RT_PROCESSIVITY_LIMIT)
print(f"Resulting cDNA: {cDNA_from_polyT}")
print(f"SUCCESSFUL CAPTURE? {'CDR3_SEQUENCE' in cDNA_from_polyT}")
print("-" * 50)

# Scenario 2: Proposed method using a TCR-specific primer
print("\nSCENARIO 2: Priming from C-Region (Proposed Method - Answer A)")
cDNA_from_C_region = simulate_rt(TCR_mRNA, C_REGION_PRIMER, RT_PROCESSIVITY_LIMIT)
print(f"Resulting cDNA: {cDNA_from_C_region}")
print(f"SUCCESSFUL CAPTURE? {'CDR3_SEQUENCE' in cDNA_from_C_region}")
print("-" * 50)

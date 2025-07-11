import re

def reverse_complement(dna_sequence):
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement[base] for base in reversed(dna_sequence))

# 1. Model the TCR mRNA sequence and experimental parameters
# These are simplified representations for demonstration purposes.
V_REGION = "GTCACACAGACAGCATATGCG" * 15  # ~315 bp
CDR3_REGION = "TGTGCCAGCAGCTTGGACAGGGGGCTGGAGCTGTTTTTT"  # 45 bp, this is the key target
J_REGION = "GGAACTCGGACCCTGAGCTGCTCGTGTACA" * 2 # ~60 bp
CONSTANT_REGION = "GCTGATGGCTACAACTTCGAC" * 25 # ~525 bp
THREE_PRIME_UTR = "GTCAGTCAGTCAGTCAGTCAG" * 15 # ~315 bp
POLYA_TAIL = "A" * 50

# The full mRNA transcript from 5' to 3'
tcr_mrna = V_REGION + CDR3_REGION + J_REGION + CONSTANT_REGION + THREE_PRIME_UTR + POLYA_TAIL
tcr_mrna_len = len(tcr_mrna)
read2_length = 225

print("--- Step-by-Step Simulation ---")
print(f"Total length of model TCR mRNA: {tcr_mrna_len} bases")
print(f"Sequencing Read 2 length: {read2_length} bases")
print("-" * 30)

# 2. Simulate the Original Failed Method
print("\nSimulation 1: Original Method (Sequencing from 3' end)")
# In a 3' capture assay, Read 2 sequences the cDNA from the polyA tail end.
# This corresponds to the last 'read2_length' bases of the mRNA sequence.
sequence_from_original_method = tcr_mrna[-read2_length:]
contains_cdr3_original = CDR3_REGION in sequence_from_original_method

print(f"Sequencing the last {read2_length} bases of the transcript.")
print(f"Does the captured sequence contain the CDR3 region? {'Yes' if contains_cdr3_original else 'No'}")
if not contains_cdr3_original:
    cdr3_pos_start = tcr_mrna.find(CDR3_REGION)
    cdr3_pos_end = cdr3_pos_start + len(CDR3_REGION)
    print(f"Reason: The CDR3 is located at bases {cdr3_pos_start}-{cdr3_pos_end} from the 5' end.")
    print(f"The 225 bp read only covers the last part of the transcript, which is the 3' UTR and PolyA tail.")
print("-" * 30)


# 3. Simulate the Proposed Solution (Option D)
print("\nSimulation 2: Proposed Method (Targeted PCR Amplification)")
print("Strategy: Use a PCR primer that binds to the TCR Constant Region.")

# Define a specific primer that binds within the Constant region.
# This would be a reverse primer for the cDNA, so it's the same sequence as the mRNA.
pcr_primer_sequence = "GCTGATGGCTACAACTTCGAC"
# Find the primer binding site on the mRNA. We'll choose a site far enough upstream.
primer_binding_site_end = tcr_mrna.find(pcr_primer_sequence, len(V_REGION) + len(CDR3_REGION) + len(J_REGION) + 50)
primer_binding_site_start = primer_binding_site_end - len(pcr_primer_sequence)

# The resulting amplicon is sequenced by Read 2 starting from this primer site and reading towards the 5' end.
start_seq = max(0, primer_binding_site_end - read2_length)
end_seq = primer_binding_site_end
sequence_from_new_method = tcr_mrna[start_seq:end_seq]
contains_cdr3_new = CDR3_REGION in sequence_from_new_method

print(f"PCR primer binds to the Constant Region at position: {primer_binding_site_end}")
print(f"Sequencing {read2_length} bases upstream from this primer site.")
print(f"Does the new captured sequence contain the CDR3 region? {'Yes' if contains_cdr3_new else 'No'}")
if contains_cdr3_new:
    print(f"Success! The amplified sequence now includes the critical V(D)J junction.")
    # Show the captured sequence context
    print("\n--- Equation of the new amplicon sequence ---")
    captured_v = V_REGION[V_REGION.find(re.search(f"({V_REGION})$", sequence_from_new_method).group(1)):]
    captured_cdr3 = CDR3_REGION
    captured_j = J_REGION
    captured_c = CONSTANT_REGION[:CONSTANT_REGION.find(sequence_from_new_method.split(J_REGION)[1])]
    
    print(f"V Region Part Length = {len(captured_v)}")
    print(f"CDR3 Region Length = {len(captured_cdr3)}")
    print(f"J Region Length = {len(captured_j)}")
    print(f"Constant Region Part Length = {len(captured_c)}")
    total_len = len(captured_v) + len(captured_cdr3) + len(captured_j) + len(captured_c)
    print(f"Total Sequenced Region = {len(captured_v)} (V) + {len(captured_cdr3)} (CDR3) + {len(captured_j)} (J) + {len(captured_c)} (C) = {total_len} bp")

<<<D>>>
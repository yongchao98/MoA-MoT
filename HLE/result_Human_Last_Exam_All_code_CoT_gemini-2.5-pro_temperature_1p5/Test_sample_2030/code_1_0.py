import textwrap

def get_reading_frames(sequence):
    """Generates codons for all 6 reading frames from a DNA sequence."""
    frames = {}
    # Forward frames
    frames['+1'] = textwrap.wrap(sequence, 3)
    frames['+2'] = textwrap.wrap(sequence[1:], 3)
    frames['+3'] = textwrap.wrap(sequence[2:], 3)
    
    # Reverse complement
    rev_seq = reverse_complement(sequence)
    
    # Reverse frames
    frames['-1'] = textwrap.wrap(rev_seq, 3)
    frames['-2'] = textwrap.wrap(rev_seq[1:], 3)
    frames['-3'] = textwrap.wrap(rev_seq[2:], 3)
    
    # Remove partial codons at the end
    for frame, codons in frames.items():
        if len(codons[-1]) < 3:
            codons.pop()
    
    return frames

def reverse_complement(dna_sequence):
    """Computes the reverse complement of a DNA sequence."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # First, complement the sequence, then reverse it
    complement = "".join(complement_map[base] for base in dna_sequence)
    rev_comp = complement[::-1]
    return rev_comp

def translate(codons):
    """Translates a list of codons into an amino acid sequence."""
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Stop', 'TAG': 'Stop',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'TGT': 'C', 'TGC': 'C', 'TGA': 'Stop', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    peptide = [codon_table.get(codon, '?') for codon in codons]
    return [aa for aa in peptide if aa != '?']

def check_snp(codon, target_codons):
    """Checks if a codon can become one of the target codons via a single SNP."""
    for target in target_codons:
        diffs = 0
        for i in range(3):
            if codon[i] != target[i]:
                diffs += 1
        if diffs == 1:
            return True
    return False

# --- Main Logic ---

# Step 1: Define original sequence and get all frame translations
original_seq_with_spaces = "5' CTT CCC CGC ACA AGT GGT 3'"
dna_sequence = "".join(original_seq_with_spaces.split()[1:-1])

print("1. Translating all 6 reading frames...")
frames_codons = get_reading_frames(dna_sequence)
frames_peptides = {frame: translate(codons) for frame, codons in frames_codons.items()}

# Create a list of all amino acids produced, ignoring 'Stop'
all_aas = []
for peptide in frames_peptides.values():
    all_aas.extend([aa for aa in peptide if aa != 'Stop'])

# Step 2: Identify the frame with two unique amino acids
print("\n2. Identifying the frame with exactly two unique amino acids...")
candidate_frames = []
for frame, peptide in frames_peptides.items():
    unique_aa_in_frame = []
    # Count how many times each AA in this frame appears in the global list
    for aa in set(peptide):
        if aa != 'Stop' and all_aas.count(aa) == 1:
            unique_aa_in_frame.append(aa)
            
    if len(unique_aa_in_frame) == 2:
        candidate_frames.append(frame)
        print(f"   - Found candidate: Frame {frame} with unique amino acids {unique_aa_in_frame}")

# Step 3: Apply SNP rules to the candidate frame(s)
print("\n3. Applying SNP rules to find the correct frame...")
target_frame_name = None
aa_properties = {
    'K': 'polar', 'W': 'nonpolar', 'N': 'polar', 'M': 'nonpolar', 'C': 'polar'
}
stop_codons = ['TAA', 'TAG', 'TGA']
cysteine_codons = ['TGT', 'TGC']

for frame_name in candidate_frames:
    peptide = frames_peptides[frame_name]
    codons = frames_codons[frame_name]
    
    unique_aas = [aa for aa in set(peptide) if aa != 'Stop' and all_aas.count(aa) == 1]
    
    polar_aa = None
    nonpolar_aa = None
    
    for aa in unique_aas:
        if aa_properties.get(aa) == 'polar':
            polar_aa = aa
        elif aa_properties.get(aa) == 'nonpolar':
            nonpolar_aa = aa

    if not (polar_aa and nonpolar_aa):
        continue

    # Get the codons for the unique amino acids
    codon_polar = codons[peptide.index(polar_aa)]
    codon_nonpolar = codons[peptide.index(nonpolar_aa)]
    
    # Check if the SNP conditions are met
    polar_to_stop = check_snp(codon_polar, stop_codons)
    nonpolar_to_cys = check_snp(codon_nonpolar, cysteine_codons)
    
    print(f"   - Checking Frame {frame_name}:")
    print(f"     - Polar AA '{polar_aa}' (codon {codon_polar}) -> Stop? {'Yes' if polar_to_stop else 'No'}")
    print(f"     - Non-polar AA '{nonpolar_aa}' (codon {codon_nonpolar}) -> Cys? {'Yes' if nonpolar_to_cys else 'No'}")

    if polar_to_stop and nonpolar_to_cys:
        target_frame_name = frame_name
        print(f"   - Frame {frame_name} is the correct frame.")
        break

# Step 4 & 5: Determine translated sequence and design the oligo
print("\n4. Designing the oligo for the modified sequence...")
if target_frame_name:
    codons_of_interest = frames_codons[target_frame_name]
    peptide_of_interest = frames_peptides[target_frame_name]
    
    # Find the position of the polar amino acid that becomes a stop codon
    polar_aa = [aa for aa in set(peptide_of_interest) if aa != 'Stop' and all_aas.count(aa) == 1 and aa_properties.get(aa) == 'polar'][0]
    stop_position = peptide_of_interest.index(polar_aa)
    
    # The sequence translated into amino acids is everything before the stop codon
    translated_codons = codons_of_interest[:stop_position]
    target_dna_for_oligo = "".join(translated_codons)
    
    # The oligo is the reverse complement
    oligo_sequence = reverse_complement(target_dna_for_oligo)
    
    print(f"   - The original peptide in frame {target_frame_name} is: {'-'.join(peptide_of_interest)}")
    print(f"   - The SNP turns the codon for '{polar_aa}' into a Stop codon.")
    print(f"   - The DNA sequence that is still translated is: {' '.join(translated_codons)} ({target_dna_for_oligo})")
    
    print("\nFinal Answer:")
    final_sequence_string = f"5' {textwrap.fill(oligo_sequence, 3).replace(' ', ' ')} 3'"
    print(f"The DNA sequence for the oligo is: {final_sequence_string}")

else:
    print("Could not find a frame that meets all the specified conditions.")

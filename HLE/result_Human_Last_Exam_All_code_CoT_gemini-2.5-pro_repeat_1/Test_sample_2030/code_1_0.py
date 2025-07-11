import sys
from io import StringIO

def solve_oligo_design():
    """
    This function solves the DNA oligo design problem by following a systematic approach:
    1. Translates the original DNA sequence in all six reading frames.
    2. Identifies the frame with two unique amino acids.
    3. Simulates the specified SNPs on the codons of these unique amino acids.
    4. Determines the DNA sequence that is translated into amino acids in the modified frame.
    5. Calculates the reverse complement of this sequence to design the binding oligo.
    """
    # Step 1: Define biological data
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

    amino_acid_properties = {
        'S': 'Polar', 'T': 'Polar', 'C': 'Polar', 'N': 'Polar', 'Q': 'Polar', 'Y': 'Polar',
        'G': 'Non-polar', 'A': 'Non-polar', 'V': 'Non-polar', 'L': 'Non-polar',
        'I': 'Non-polar', 'P': 'Non-polar', 'F': 'Non-polar', 'M': 'Non-polar', 'W': 'Non-polar',
    }

    cysteine_codons = {'TGT', 'TGC'}
    stop_codons = {'TAA', 'TAG', 'TGA'}
    original_seq = "CTTCCCGCACAAGTGGT"

    def reverse_complement(dna_seq):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return "".join(complement.get(base, base) for base in reversed(dna_seq))

    def translate_frame(dna_seq, frame_start):
        peptide = ""
        codons = []
        for i in range(frame_start, len(dna_seq) - 2, 3):
            codon = dna_seq[i:i+3]
            codons.append(codon)
            peptide += codon_table.get(codon, '')
        return peptide, codons

    # Step 2: Six-frame translation
    translations = {}
    reverse_seq = reverse_complement(original_seq)
    for i in range(3):
        p, c = translate_frame(original_seq, i)
        translations[f"Frame +{i+1}"] = {'peptide': p, 'codons': c}
    for i in range(3):
        p, c = translate_frame(reverse_seq, i)
        translations[f"Frame -{i+1}"] = {'peptide': p, 'codons': c}

    # Step 3: Identify the target frame
    target_frame_name = None
    unique_aas_in_frame = []
    for frame_name, data in translations.items():
        other_peptides_set = set()
        for other_frame, other_data in translations.items():
            if frame_name != other_frame:
                other_peptides_set.update(list(other_data['peptide']))
        
        current_unique = [aa for aa in set(data['peptide']) if aa not in other_peptides_set]
        if len(current_unique) == 2:
            target_frame_name = frame_name
            unique_aas_in_frame = current_unique
            break
            
    print(f"Found target frame: {target_frame_name}")
    print(f"Unique amino acids in this frame: {', '.join(unique_aas_in_frame)}")
    
    # Step 4: Simulate SNPs
    target_data = translations[target_frame_name]
    peptide_seq = target_data['peptide']
    codon_seq = target_data['codons']

    polar_aa = None
    nonpolar_aa = None
    if amino_acid_properties.get(unique_aas_in_frame[0]) == 'Polar':
        polar_aa, nonpolar_aa = unique_aas_in_frame[0], unique_aas_in_frame[1]
    else:
        polar_aa, nonpolar_aa = unique_aas_in_frame[1], unique_aas_in_frame[0]

    polar_aa_codon = codon_seq[peptide_seq.find(polar_aa)]
    nonpolar_aa_codon = codon_seq[peptide_seq.find(nonpolar_aa)]

    def find_snp(from_codon, to_codons):
        for to_codon in to_codons:
            if sum(1 for i in range(3) if from_codon[i] != to_codon[i]) == 1:
                return to_codon
        return None

    modified_polar_codon = find_snp(polar_aa_codon, stop_codons)
    modified_nonpolar_codon = find_snp(nonpolar_aa_codon, cysteine_codons)

    print(f"\nSNP analysis:")
    print(f"1. Polar amino acid '{polar_aa}' (codon {polar_aa_codon}) is changed to a Stop codon: {modified_polar_codon}")
    print(f"2. Non-polar amino acid '{nonpolar_aa}' (codon {nonpolar_aa_codon}) is changed to Cysteine: {modified_nonpolar_codon}")

    # Step 5: Construct modified sequence and design oligo
    modified_codons = [c.replace(polar_aa_codon, modified_polar_codon).replace(nonpolar_aa_codon, modified_nonpolar_codon) for c in codon_seq]
    
    translated_codons = []
    for codon in modified_codons:
        if codon in stop_codons:
            break
        translated_codons.append(codon)
    
    modified_translated_dna = "".join(translated_codons)
    oligo_sequence = reverse_complement(modified_translated_dna)

    # Step 6: Output the result
    print(f"\nModified DNA sequence translated into amino acids: 5' {modified_translated_dna} 3'")
    print(f"The oligo must bind to this sequence.")
    print(f"\nFinal oligo sequence (reverse complement):")
    final_answer = f"5' {oligo_sequence} 3'"
    print(final_answer)
    
    # Capture the final sequence for the specific format requirement
    # Use a string buffer to capture print output if needed, but for this problem, direct return is cleaner.
    return oligo_sequence

# Execute the solver function
final_oligo_seq = solve_oligo_design()
print(f"\n<<<5' {final_oligo_seq} 3'>>>")

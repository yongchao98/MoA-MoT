import sys
import io

def check_correctness():
    """
    This function checks the correctness of the provided LLM answer by analyzing the genetic sequences.
    The goal is to find the mutation most likely to create a non-functional protein, which would
    eliminate the anti-nutritional compound. The most effective mutations for this purpose are
    nonsense mutations (creating a premature stop codon) or frameshift mutations.
    """

    # Standard DNA codon table mapping codons to one-letter amino acid codes. '*' represents a stop codon.
    CODON_TABLE = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    def translate(seq):
        """Translates a DNA sequence into an amino acid sequence until a stop codon is found."""
        protein = []
        # Ensure the sequence length is a multiple of 3 for codon processing
        seq_len = len(seq) - (len(seq) % 3)
        for i in range(0, seq_len, 3):
            codon = seq[i:i+3]
            aa = CODON_TABLE.get(codon, 'X')  # 'X' for unknown codons
            protein.append(aa)
            if aa == '*': # Stop translation after a stop codon
                break
        return "".join(protein)

    def get_mutation_type_and_disruptiveness(original_seq, mutant_seq):
        """
        Analyzes a mutation and assigns a numerical score for its disruptiveness.
        Higher score means more disruptive.
        3: Highly disruptive (Nonsense, Frameshift)
        2: Disruptive (In-frame indel)
        1: Potentially disruptive (Missense)
        0: Not disruptive (Silent)
        """
        # Check for frameshift (indel not a multiple of 3)
        len_diff = len(mutant_seq) - len(original_seq)
        if len_diff != 0 and len_diff % 3 != 0:
            return "Frameshift", 3

        original_protein = translate(original_seq)
        mutant_protein = translate(mutant_seq)

        # Check for nonsense mutation (premature stop)
        if '*' in mutant_protein and mutant_protein.find('*') < len(original_protein) - 1:
            return "Nonsense", 3

        # Check for in-frame indel (insertion/deletion of whole codons)
        if len(original_protein) != len(mutant_protein):
            return "In-frame indel", 2

        # Check for missense mutation
        if original_protein != mutant_protein:
            return "Missense", 1

        # Check for silent mutation
        if original_seq != mutant_seq:
            return "Silent", 0

        return "No change", 0

    # Define sequences from the question
    intact_gene = "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT"
    mutants = {
        "Mutant 1": "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 2": "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 3": "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT",
        "Mutant 4": "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT",
    }
    
    # The LLM's answer is D, which corresponds to Mutant 2
    llm_answer_choice = "Mutant 2"

    # Analyze each mutant to find the most disruptive one
    disruption_scores = {}
    mutation_types = {}
    for name, seq in mutants.items():
        mutation_type, score = get_mutation_type_and_disruptiveness(intact_gene, seq)
        disruption_scores[name] = score
        mutation_types[name] = mutation_type

    # Find the mutant with the highest disruption score
    most_disruptive_mutant = max(disruption_scores, key=disruption_scores.get)

    # Check if the LLM's choice matches the most disruptive mutant
    if llm_answer_choice == most_disruptive_mutant:
        # Verify the specific reason for the choice is correct
        m2_type = mutation_types["Mutant 2"]
        if m2_type == "Nonsense":
            return "Correct"
        else:
            return f"Incorrect. The answer D (Mutant 2) is chosen, but the reason is flawed. The mutation is a {m2_type}, not a nonsense mutation as claimed."
    else:
        m_type = mutation_types[most_disruptive_mutant]
        return f"Incorrect. The provided answer is D (Mutant 2), but the most disruptive mutation is {most_disruptive_mutant}, which causes a '{m_type}' mutation. A more disruptive mutation is more likely to eliminate the compound."

# Execute the check and return the result
result = check_correctness()
print(result)
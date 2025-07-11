import collections

def predict_oligomeric_states(sequences):
    """
    Predicts the oligomeric state of coiled-coil sequences based on dominant
    core 'a' and 'd' residues.
    """
    
    # Rules based on seminal protein design literature for (a, d) pairs.
    # I=Isoleucine, A=Alanine, L=Leucine, Q=Glutamine, T=Threonine
    prediction_rules = {
        ('I', 'A'): 2,  # Dimer
        ('I', 'L'): 3,  # Trimer
        ('I', 'I'): 4,  # Tetramer
        ('I', 'Q'): 5,  # Pentamer
        ('I', 'T'): 7,  # Heptamer
    }
    
    oligomeric_states = []
    
    for seq in sequences:
        a_residues = []
        d_residues = []
        
        # Extract 'a' position residues (indices 1, 8, 15, 22...)
        for i in range(1, len(seq), 7):
            a_residues.append(seq[i])
            
        # Extract 'd' position residues (indices 4, 11, 18, 25...)
        for i in range(4, len(seq), 7):
            d_residues.append(seq[i])
            
        # Determine the dominant (most common) residue at each position
        # This handles minor mutations within the sequence
        dominant_a = collections.Counter(a_residues).most_common(1)[0][0]
        dominant_d = collections.Counter(d_residues).most_common(1)[0][0]
        
        core_pair = (dominant_a, dominant_d)
        
        # Look up the oligomeric state from the rules
        state = prediction_rules.get(core_pair, "Unknown")
        oligomeric_states.append(state)
        
    # The prompt asks to output each number in the final equation.
    # We will print the sequence of predicted states.
    # Sequence 1: ('I', 'A') -> 2
    # Sequence 2: ('I', 'L') -> 3
    # Sequence 3: ('I', 'I') -> 4
    # Sequence 4: ('I', 'Q') -> 5
    # Sequence 5: ('I', 'T') -> 7
    # Final Result: 2, 3, 4, 5, 7
    
    print("The predicted oligomeric states are:")
    # Print each number, separated by commas
    print(*oligomeric_states, sep=", ")

# The list of coiled-coil sequences to analyze
protein_sequences = [
    "EIAQALKEIAKALKEIAWALKEIAQALK",
    "EIAALKQEIAALKKENAALKQEIAALKQ",
    "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
    "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
    "EIAQTLKEIAKTLKEIAWTLKEIAQTLK"
]

predict_oligomeric_states(protein_sequences)

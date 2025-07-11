import collections

def predict_coiled_coil_state(sequences):
    """
    Predicts the oligomeric state of coiled-coil sequences based on the
    residues at the 'a' and 'd' positions of the heptad repeat.
    """
    
    predictions = []
    
    print("Analyzing coiled-coil sequences...")
    print("-" * 35)

    for seq in sequences:
        # Assuming a 'gabcdef' register where the first Isoleucine is at 'a'.
        # 'a' positions are 1, 8, 15, 22...
        # 'd' positions are 4, 11, 18, 25...
        a_residues = [seq[i] for i in range(1, len(seq), 7)]
        d_residues = [seq[i] for i in range(4, len(seq), 7)]
        
        # Determine the most common 'a' and 'd' residues to define the core pair
        a_res = collections.Counter(a_residues).most_common(1)[0][0]
        d_res = collections.Counter(d_residues).most_common(1)[0][0]
        
        state = 'Unknown'
        
        # Prediction rules based on common (a, d) pairings
        if a_res == 'I' and d_res == 'A':
            # IaAd packing is characteristic of trimers.
            state = 3
        elif a_res == 'I' and d_res == 'L':
            # IaLd is the canonical pairing for dimers. The one 'N' in sequence 2
            # is not enough to override the strong dimeric signal.
            state = 2
        elif a_res == 'I' and d_res == 'I':
            # IaId packing is characteristic of trimers.
            state = 3
        elif a_res == 'I' and d_res == 'Q':
            # IaQd, with a polar residue, can be accommodated in a trimeric core.
            state = 3
        elif a_res == 'I' and d_res == 'T':
            # IaTd has beta-branched side chains at both positions, which can
            # favor a tetrameric arrangement.
            state = 4
            
        predictions.append(state)
        print(f"Sequence: {seq}")
        print(f"Core pair (a,d): ({a_res}, {d_res})")
        print(f"Predicted oligomeric state: {state}\n")
        
    # Print the final combined result
    final_result_str = ", ".join(map(str, predictions))
    print("-" * 35)
    print(f"The final predicted oligomeric states are: {final_result_str}")


# The provided coiled-coil sequences
protein_sequences = [
    "EIAQALKEIAKALKEIAWALKEIAQALK",
    "EIAALKQEIAALKKENAALKQEIAALKQ",
    "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
    "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
    "EIAQTLKEIAKTLKEIAWTLKEIAQTLK"
]

predict_coiled_coil_state(protein_sequences)
<<<E>>>
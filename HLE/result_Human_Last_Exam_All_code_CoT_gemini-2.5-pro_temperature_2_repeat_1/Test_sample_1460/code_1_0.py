def reduce_braid_word(braid_word, relevant_strands):
    """
    Reduces a braid word from a larger braid group to a smaller one
    based on a subset of relevant strands.
    
    Args:
        braid_word (str): The braid word in B_n, e.g., "s1^2 s2^2 s3 s4^-1".
        relevant_strands (list): A sorted list of the 1-indexed strands to focus on.
        
    Returns:
        A tuple containing the new braid word in B_k and the size k of the new group.
    """
    k = len(relevant_strands)
    strand_map = {original_strand: new_strand for new_strand, original_strand in enumerate(relevant_strands, 1)}
    
    new_braid_terms = []
    
    terms = braid_word.split()
    for term in terms:
        # Parse term like "s_i^p"
        if '^' in term:
            base, power = term.split('^')
            power = int(power)
        else:
            base, power = term, 1
            
        original_index = int(base[1:])
        
        # Check if this generator acts on the relevant strands
        # sigma_i acts on strands i and i+1
        if original_index in strand_map and original_index + 1 in strand_map:
            new_index = strand_map[original_index]
            if power == 1:
                new_braid_terms.append(f"s{new_index}")
            elif power == -1:
                new_braid_terms.append(f"s{new_index}^-1")
            else:
                new_braid_terms.append(f"s{new_index}^{power}")

    return " ".join(new_braid_terms), k

# The original braid word in B_5
original_braid_word = "s1^2 s2^2 s3 s4^-1"
# The strands forming the knot we want to identify
component_strands = [3, 4, 5]

# The part of the word that braids the component strands amongst themselves
# We ignore s2^2 as it only creates linking with another component
word_for_knot_type = "s3 s4^-1"

# Reduce the problem to a smaller braid group
final_braid, new_n = reduce_braid_word(word_for_knot_type, component_strands)

print(f"The original braid is: β = σ₁²σ₂²σ₃σ₄⁻¹ in B₅")
print(f"The connected component in question is formed by strands {component_strands}.")
print(f"The braiding that determines the knot type of this component is σ₃σ₄⁻¹.")
print(f"By re-indexing strands {component_strands} to be {{1, 2, ..., {new_n}}}, this corresponds to a braid in B{new_n}.")
print(f"The equivalent braid in B{new_n} is: {final_braid.replace('s', 'σ').replace('^', '')}")
print("The closure of the braid σ₁σ₂⁻¹ is known in knot theory to be the Trefoil knot.")

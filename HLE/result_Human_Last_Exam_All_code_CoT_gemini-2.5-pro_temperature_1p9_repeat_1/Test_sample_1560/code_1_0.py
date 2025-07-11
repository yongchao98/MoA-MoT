import math
from collections import defaultdict

def solve_quantum_state():
    """
    This function calculates the final state of the entangled quantum system
    after applying gate Q to the first ququint.
    """

    # The action of gate Q is defined as:
    # Q |i> = 1/sqrt(2) * ( |j> + |k> )
    # We can represent this with a dictionary mapping i to [j, k].
    q_map = {
        0: [1, 2],
        1: [0, 3],
        2: [1, 4],
        3: [2, 0],
        4: [3, 2],
    }

    # We start with the entangled state: |Ψ_AB⟩ = 1/√5 * Σ_i |i⟩_A ⊗ |i⟩_B
    # We apply Q to ququint A: |Φ_AB⟩ = (Q ⊗ I) |Ψ_AB⟩
    # |Φ_AB⟩ = 1/√5 * Σ_i (Q|i⟩_A) ⊗ |i⟩_B
    # Substituting the definition of Q:
    # |Φ_AB⟩ = 1/√5 * Σ_i [ 1/√2 * (|j>_A + |k>_A) ] ⊗ |i⟩_B
    # |Φ_AB⟩ = 1/√10 * Σ_i ( |j>_A|i>_B + |k>_A|i>_B ) where [j, k] = q_map[i]

    # Let's collect terms by the state of ququint A (|k>_A).
    # The dictionary will map an index k to a list of indices l for state |k>A ⊗ |l>B.
    grouped_terms = defaultdict(list)
    
    for i in range(5):
        # The initial state term is |i>A ⊗ |i>B
        original_ket_B_idx = i
        
        # Applying Q to |i>A results in two terms
        new_ket_A_idx_1 = q_map[i][0]
        new_ket_A_idx_2 = q_map[i][1]
        
        # Add the resulting B-kets to the lists for the corresponding A-kets
        grouped_terms[new_ket_A_idx_1].append(original_ket_B_idx)
        grouped_terms[new_ket_A_idx_2].append(original_ket_B_idx)

    # Now, build the final equation string for printing.
    # The normalization constant is 1/sqrt(10).
    final_equation = f"|final_state> = (1/sqrt({10})) * [\n"
    
    term_lines = []
    # Sort by the A-ket index for a canonical representation.
    for k in sorted(grouped_terms.keys()):
        b_kets_indices = grouped_terms[k]
        
        a_part_str = f"|{k}>A"
        
        # Format the part for ququint B: (|l1>B + |l2>B + ...)
        # Sort B-ket indices for a canonical representation.
        b_kets_sorted_indices = sorted(b_kets_indices)
        if len(b_kets_sorted_indices) == 1:
            # e.g., |2>B
            b_part_str = f"|{b_kets_sorted_indices[0]}>B"
        else:
            # e.g., (|0>B + |3>B)
            b_part_str = "(" + " + ".join([f"|{l}>B" for l in b_kets_sorted_indices]) + ")"

        # Combine the parts for this term using a tensor product symbol
        term_line = f"  |{k}>A ⊗ {b_part_str}"
        term_lines.append(term_line)
        
    final_equation += " +\n".join(term_lines)
    final_equation += "\n]"

    print(final_equation)
    # The final equation is the content for the answer
    return final_equation

# Execute the function to get the final state
final_answer = solve_quantum_state()
#<<<|final_state> = (1/sqrt(10)) * [
#  |0>A ⊗ (|1>B + |3>B) +
#  |1>A ⊗ (|0>B + |2>B) +
#  |2>A ⊗ (|0>B + |3>B + |4>B) +
#  |3>A ⊗ (|1>B + |4>B) +
#  |4>A ⊗ |2>B
#]>>>
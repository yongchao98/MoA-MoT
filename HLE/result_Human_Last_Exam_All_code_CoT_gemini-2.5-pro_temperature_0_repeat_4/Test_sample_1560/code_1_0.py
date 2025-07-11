import math
from fractions import Fraction
from collections import defaultdict

def solve_ququint_problem():
    """
    Calculates and prints the possible final states of a two-ququint system
    after a gate operation and measurement.
    """
    # 1. Define the action of gate Q on the basis states |i>
    # Q|i> = 1/sqrt(2) * sum of other basis states
    # We represent this as a dictionary mapping input basis index to a list of output basis indices.
    q_action = {
        0: [1, 2],
        1: [0, 3],
        2: [1, 4],
        3: [2, 0],
        4: [3, 2],
    }
    
    # The coefficient for each term in the Q output is 1/sqrt(2)
    q_coeff = 1 / math.sqrt(2)

    # 2. Define the initial entangled state
    # |Psi> = 1/sqrt(5) * sum_{i=0 to 4} |i>_A |i>_B
    # We can represent this as a list of tuples (coeff, basis_A, basis_B)
    initial_coeff = 1 / math.sqrt(5)
    
    # 3. Apply the gate Q to ququint A
    # The new state will have terms with coefficient initial_coeff * q_coeff = 1/sqrt(10)
    post_q_state = []
    for i in range(5):
        # Get the result of Q|i>
        output_bases_A = q_action[i]
        for basis_A in output_bases_A:
            # The state is |basis_A>_A |i>_B with a coefficient of 1/sqrt(10)
            post_q_state.append({'coeff': initial_coeff * q_coeff, 'A': basis_A, 'B': i})

    # 4. Regroup terms by the state of ququint A to prepare for measurement
    regrouped_by_A = defaultdict(list)
    for term in post_q_state:
        regrouped_by_A[term['A']].append(term['B'])

    print("The gate Q is applied to the first ququint of the entangled pair, and then it is measured.")
    print("The final state of the system depends on the measurement outcome. Here are all possibilities:\n")

    # 5. Calculate and print the outcome for each possible measurement of A
    for i in sorted(regrouped_by_A.keys()):
        b_states = regrouped_by_A[i]
        num_terms = len(b_states)
        
        # Probability of measuring |i>_A is (num_terms * |1/sqrt(10)|^2)
        prob = Fraction(num_terms, 10)
        
        # The state of ququint B collapses to a superposition of the corresponding states.
        # The normalization coefficient for this new state is 1/sqrt(num_terms).
        norm_coeff_B = 1 / math.sqrt(num_terms)
        
        # Build the string representation of the final state
        b_state_str_parts = []
        for basis_B in sorted(b_states):
            b_state_str_parts.append(f"|{basis_B}>_B")
        
        b_superposition_str = " + ".join(b_state_str_parts)

        if num_terms == 1:
            # If there's only one term, the state is not a superposition
            final_b_state_str = b_superposition_str
        else:
            # For superpositions, show the coefficient
            norm_val = int(round(1 / norm_coeff_B**2))
            final_b_state_str = f"(1/sqrt({norm_val})) * ({b_superposition_str})"

        print(f"Outcome {i}: Ququint A is measured as |{i}>_A")
        print(f"  - Probability: {prob.numerator}/{prob.denominator}")
        # Final equation with all numbers
        print(f"  - Final State: |{i}>_A ⊗ {final_b_state_str}\n")

solve_ququint_problem()
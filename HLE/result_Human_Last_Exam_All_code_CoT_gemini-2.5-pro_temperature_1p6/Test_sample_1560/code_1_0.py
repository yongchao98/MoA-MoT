import math
from collections import defaultdict

def solve_quantum_problem():
    """
    This script solves the ququint problem by applying the gate Q to the first 
    ququint of an entangled pair and determines the possible final states 
    of the second ququint after a measurement on the first.
    """

    # 1. Define the action of gate Q on basis states |i>.
    # The definition Q|i> = 1/sqrt(2) * (sum of two states) is used.
    # We represent the output as a list of basis state indices.
    q_action = {
        0: [1, 2],  # Q|0> -> 1/sqrt(2) * (|1> + |2>)
        1: [0, 3],  # Q|1> -> 1/sqrt(2) * (|0> + |3>)
        2: [1, 4],  # Q|2> -> 1/sqrt(2) * (|1> + |4>)
        3: [2, 0],  # Q|3> -> 1/sqrt(2) * (|2> + |0>)
        4: [3, 2],  # Q|4> -> 1/sqrt(2) * (|3> + |2>)
    }
    
    # The normalization factor for the amplitude of each basis tensor |A,B>
    # comes from 1/sqrt(5) (initial state) and 1/sqrt(2) (gate Q),
    # resulting in a total factor of 1/sqrt(10). The squared norm is 1/10.
    norm_factor_sq_per_term = 10

    # 2. Start with the initial entangled state.
    # |Psi_ent> = 1/sqrt(5) * sum_{i=0 to 4} |i>_A |i>_B
    # This is represented as a list of (state_A, state_B) pairs.
    initial_pairs = [(i, i) for i in range(5)]

    # 3. Apply Q to the first ququint (A).
    # The operation is (Q x I) |i>_A |i>_B = (Q|i>_A) |i>_B.
    # We generate a list of all resulting |A,B> tensor product states.
    final_state_terms = []
    for state_a, state_b in initial_pairs:
        result_a_states = q_action[state_a]
        for new_state_a in result_a_states:
            final_state_terms.append((new_state_a, state_b))
    
    # 4. Group terms by the state of ququint A to analyze measurement outcomes.
    measurement_map = defaultdict(list)
    for state_a, state_b in final_state_terms:
        measurement_map[state_a].append(state_b)
        
    print("When gate Q is applied to the first ququint of the entangled pair, and a measurement is performed on it, the system collapses.")
    print("The final state of the second ququint depends on the probabilistic measurement outcome of the first.\n")
    print("The possible outcomes are:\n")

    # 5. For each possible measurement outcome on A, calculate the probability
    # and the resulting normalized state of ququint B.
    final_answer_candidate = ""
    for outcome_a in sorted(measurement_map.keys()):
        # These are the basis states of B associated with outcome_a
        b_states = measurement_map[outcome_a]
        
        # Probability is the number of terms for this outcome divided by the total number of terms.
        num_terms = len(b_states)
        probability_numerator = num_terms
        probability_denominator = norm_factor_sq_per_term
        
        # To find the normalized state of B, we count coefficients.
        counts = defaultdict(int)
        for state in b_states:
            counts[state] += 1
        
        # The normalization factor for B's state vector is sqrt(sum of |coeffs|^2),
        # which is sqrt(num_terms) since all unnormalized coefficients are 1.
        norm_b_sq = num_terms
        
        # Build the string representation for the state of B.
        state_b_parts = []
        for state in sorted(counts.keys()):
            coeff = counts[state]
            if norm_b_sq == 1:
                 # If norm is 1, coefficient is just 1.
                 state_b_parts.append(f"|{state}>_B")
            else:
                 # The prompt asks for each number in the equation.
                 state_b_parts.append(f"({coeff}/sqrt({norm_b_sq}))|{state}>_B")
        
        state_b_str = " + ".join(state_b_parts)

        # Print the detailed result for this measurement outcome.
        print(f"If ququint A is measured as |{outcome_a}> (with probability {probability_numerator}/{probability_denominator}):")
        print(f"  The state of ququint B collapses to: {state_b_str}")
        print("-" * 30)

        # Storing the most complex state as a candidate for the final answer block.
        if outcome_a == 2:
            # Let's rebuild the string for the specific required format without the "_B" part.
            final_eq_parts = []
            for state in sorted(counts.keys()):
                 coeff = counts[state]
                 final_eq_parts.append(f"({coeff}/sqrt({norm_b_sq}))|{state}>")
            final_answer_candidate = " + ".join(final_eq_parts)

    # Return the candidate for the final answer block.
    return final_answer_candidate


final_state_equation = solve_quantum_problem()
print(f"\n<<<As an example, the final state equation for ququint B if ququint A is measured as |2> is: {final_state_equation}>>>")

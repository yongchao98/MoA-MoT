import math

def describe_measurement_outcomes():
    """
    Analyzes and describes the possible outcomes of measuring one ququint
    from a two-ququint entangled pair.
    """
    num_states = 5
    
    # The initial entangled state is:
    # 1/sqrt(5) * ( |0>_A|0>_B + |1>_A|1>_B + ... + |4>_A|4>_B )
    # The coefficient for each basis pair is 1/sqrt(5).
    coefficient = 1 / math.sqrt(num_states)
    
    # The probability of measuring any specific state |k> on ququint A is
    # the square of the magnitude of its coefficient.
    probability_val = coefficient**2
    
    print("The initial entangled state of the two-ququint system is:")
    initial_state_str = " + ".join([f"|{i}>_A ⊗ |{i}>_B" for i in range(num_states)])
    print(f"Ψ = (1/√{num_states}) * ({initial_state_str})\n")
    
    print("When ququint A is measured in the computational basis, there are 5 possible outcomes.")
    print("The final state of the system depends on the measurement result.")
    print("-" * 70)
    
    # Loop through each possible measurement outcome
    for k in range(num_states):
        
        # The equation for the final state after collapse is: Final State = 1 * |k>_A ⊗ |k>_B
        # The numbers in the equation are the coefficient (1) and the state labels (k).
        final_state_equation = f"1 * |{k}>_A ⊗ |{k}>_B"
        
        print(f"Possibility {k + 1}:")
        print(f"  - If measurement on Ququint A yields outcome: |{k}>")
        print(f"  - The probability of this outcome is: 1/{num_states} (or {probability_val:.1f})")
        print(f"  - The system collapses to the final state described by the equation:")
        print(f"    Final State = {final_state_equation}")
        print("-" * 70)

# Run the analysis
describe_measurement_outcomes()

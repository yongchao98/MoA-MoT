import textwrap

def explain_recombination_frequency():
    """
    Explains why the highest recombination frequency occurs at the earliest
    transferred gene in an interrupted mating experiment.
    """
    
    # Gene order provided in the problem. 'thr' is transferred first.
    gene_order = ['thr', 'azi', 'gal']

    # Principle of Interrupted Mating
    principle = (
        "In an E. coli interrupted mating experiment, the donor chromosome is transferred linearly "
        "to the recipient cell over time. The process can be interrupted spontaneously at any moment. "
        "Therefore, genes that are closer to the origin of transfer and enter the recipient cell "
        "earlier have a higher probability of being successfully transferred than genes that are "
        "transferred later."
    )

    # Applying the principle to the given gene order
    analysis = (
        f"The given gene transfer order is {gene_order[0]} -> {gene_order[1]} -> {gene_order[2]}. "
        f"The '{gene_order[0]}' marker is the first to enter the recipient cell. "
        f"The '{gene_order[1]}' marker enters next, followed by the '{gene_order[2]}' marker. "
        "Because of the increasing chance of interruption over time, the frequency of transfer "
        f"will be highest for '{gene_order[0]}', lower for '{gene_order[1]}', and lowest for '{gene_order[2]}'."
    )
    
    # Conclusion
    conclusion = (
        "The frequency of recombinants in the recipient population directly reflects the "
        "frequency of gene transfer. Consequently, the highest frequency of recombinants will be "
        f"observed for the earliest marker transferred, which is '{gene_order[0]}'."
    )

    # Print the explanation
    print("### Analysis of Recombination Frequency ###\n")
    print(textwrap.fill(principle, width=80))
    print("\n" + "="*40 + "\n")
    print(textwrap.fill(analysis, width=80))
    print("\n" + "="*40 + "\n")
    print("Conclusion:")
    print(textwrap.fill(conclusion, width=80))
    print("\nMatching this conclusion to the options, the highest frequency of recombinants should be expected 'Immediately after thr+'.")

# Execute the explanation function
explain_recombination_frequency()
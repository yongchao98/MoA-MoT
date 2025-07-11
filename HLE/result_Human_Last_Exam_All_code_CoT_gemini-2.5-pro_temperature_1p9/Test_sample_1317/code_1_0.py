import textwrap

def explain_recombination_frequency():
    """
    Explains the reasoning behind recombination frequency in an interrupted mating experiment.
    """
    explanation_title = "Analysis of Recombination Frequency in Interrupted Mating"
    principle_header = "1. Principle of Interrupted Mating:"
    principle_text = (
        "In bacterial conjugation, an Hfr (High-Frequency Recombination) donor cell transfers "
        "its chromosome to a recipient cell in a linear, sequential manner starting from an "
        "origin of transfer (OriT). The longer the cells are allowed to mate, the more of the "
        "chromosome is transferred. The frequency of recombination for a specific gene is "
        "directly related to how early it is transferred. Earlier genes have a higher probability "
        "of being successfully transferred and integrated before the mating process is randomly interrupted."
    )

    application_header = "2. Application to the Given Scenario:"
    application_text_template = (
        "The problem states the known gene order is '{gene_order}' and the observed transfer order "
        "is that '{first_gene}' is transferred first. This means the origin of transfer (OriT) is "
        "located just before the '{first_gene}' gene on the chromosome. The transfer proceeds from "
        "there, so the order of transfer is: OriT -> {gene_order_arrow}."
    )
    gene_order = "thr-azi-gal"
    first_gene = "thr+"
    gene_order_arrow = "thr -> azi -> gal"
    application_text = application_text_template.format(
        gene_order=gene_order,
        first_gene=first_gene,
        gene_order_arrow=gene_order_arrow
    )


    conclusion_header = "3. Conclusion:"
    conclusion_text = (
        "Because recombination frequency is highest for the genes transferred earliest, the frequency "
        "will be greatest for 'thr', lower for 'azy', and lowest for 'gal'. Therefore, the genetic location "
        "with the highest frequency of recombinants is the one corresponding to the first gene transferred, 'thr'."
        "Among the choices, 'Immediately after thr+' best describes this region, as it is the first marked location to enter the recipient cell following the start of transfer."
    )

    print(f"{explanation_title}\n{'-'*len(explanation_title)}")

    print(f"\n{principle_header}")
    print(textwrap.fill(principle_text, width=80))

    print(f"\n{application_header}")
    print(textwrap.fill(application_text, width=80))

    print(f"\n{conclusion_header}")
    print(textwrap.fill(conclusion_text, width=80))

explain_recombination_frequency()
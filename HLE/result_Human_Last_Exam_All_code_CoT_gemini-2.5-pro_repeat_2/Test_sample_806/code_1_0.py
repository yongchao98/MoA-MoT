import textwrap

def explain_reaction():
    """
    This function prints an explanation of the biochemical reaction that led to the patient's skin blisters.
    """
    print("The clinical scenario describes a cascade of events likely resulting from a severe cutaneous adverse reaction (SCAR) like Stevens-Johnson Syndrome (SJS), triggered by the third medication.")
    print("\nThe biochemical reaction that initiated this process is a specific type of immune system activation detailed below:")

    steps = {
        "1. Drug-Protein Binding": "The process begins with the non-covalent binding of the drug molecule (or one of its metabolites) directly to a specific Human Leukocyte Antigen (HLA) protein. This protein is located on the surface of the patient's own cells.",
        "2. Altered 'Self' Presentation": "This binding alters the shape of the HLA molecule. The immune system, specifically T-cells, now perceives this altered HLA-drug complex as a foreign or dangerous signal.",
        "3. T-Cell Activation": "The T-cell receptor (TCR) on a cytotoxic T-lymphocyte directly recognizes and binds to this unique drug-HLA complex. This interaction bypasses the conventional antigen-processing pathway and immediately activates the T-cell.",
        "4. Cytotoxic Molecule Release": "Upon activation, the T-cell releases large quantities of cytotoxic agents. The key protein released is Granulysin.",
        "5. Widespread Cell Death": "Granulysin induces widespread apoptosis (programmed cell death) in keratinocytes (skin cells). This mass cell death causes the epidermis to detach from the dermis, leading to the formation of the blisters and skin sloughing characteristic of SJS."
    }

    for title, description in steps.items():
        print(f"\n{title}:")
        # Wrap the text for better readability
        wrapped_text = textwrap.fill(description, width=80)
        print(wrapped_text)

    print("\n" + "-"*40)
    summary = "In summary, the single biochemical reaction that initiated the entire process was the direct binding of the drug to the HLA protein, which triggered the fatal immune response against the patient's own skin cells."
    print(textwrap.fill(summary, width=80))

explain_reaction()
import textwrap

def identify_compound_1():
    """
    Identifies Compound 1 based on the reaction of geraniol with
    O-(p-tolyl) chlorothionoformate and explains the NMR data.
    """
    # Reactant Information
    geraniol = "Geraniol [(2E)-3,7-dimethylocta-2,6-dien-1-ol]"
    reagent = "O-(p-tolyl) chlorothionoformate"
    
    # Product Information
    compound_1_name = "O-((2E)-3,7-dimethylocta-2,6-dien-1-yl) O'-(p-tolyl) carbonothioate"
    compound_1_common_name = "O-geranyl O'-p-tolyl thionocarbonate"
    compound_1_formula = "C18H24O2S"
    
    # NMR Data from the problem
    geraniol_peak = "5.32-5.37 ppm (multiplet, 1H)"
    compound_1_peak = "5.97 ppm (doublet of doublets, 1H)"
    
    # Explanation
    explanation = f"""
    The reaction between the alcohol group of geraniol and {reagent} is a nucleophilic acyl substitution, which forms a thionocarbonate ester.

    Compound 1 is: {compound_1_name}
    Also known as: {compound_1_common_name}
    Chemical Formula: {compound_1_formula}

    --- NMR Rationale ---

    1.  The peak in geraniol at {geraniol_peak} corresponds to the vinylic proton at the C2 position.

    2.  Downfield Shift: This proton's signal shifts downfield to {compound_1_peak} in Compound 1. This is because the newly attached thionocarbonate group [-O-C(=S)-O-p-tolyl] has a C=S bond that exerts a strong deshielding (anisotropic) effect on the nearby C2 proton.

    3.  Splitting Pattern Change: The splitting changes from a multiplet to a doublet of doublets. This is because the bulky new group restricts free rotation around the C1-C2 single bond. This makes the two protons on the C1 carbon (-CH2-) magnetically inequivalent. The C2 proton now 'sees' two different neighboring protons, splitting its signal into a doublet of doublets.
    """
    
    print(textwrap.dedent(explanation).strip())

if __name__ == '__main__':
    identify_compound_1()
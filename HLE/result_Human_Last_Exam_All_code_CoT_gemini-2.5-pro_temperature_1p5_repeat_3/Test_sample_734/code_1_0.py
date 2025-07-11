import textwrap

def theological_answer():
    """
    This function provides the answer to the user's question about
    the Christian theologian who first developed terminology for
    changes in the laws of nature.
    """
    theologian = "Augustine of Hippo"

    explanation = (
        "Augustine of Hippo (354–430 AD) is the theologian who first developed specialized terminology "
        "to describe apparent changes in the laws of nature. He proposed the concept of 'seminal reasons' "
        "or 'causal reasons' (Latin: rationes seminales or rationes causales), which he argued God "
        "had placed within creation from the very beginning.\n\n"
        "This framework allowed him to explain:\n"
        "1.  The Virgin Birth: Not as a violation of a natural order God had to break, but as the "
        "activation of a pre-ordained potentiality that existed outside the 'ordinary course of nature' "
        "(naturae cursus).\n"
        "2.  The Fall: Changes in the natural world (e.g., the earth bringing forth 'thorns and thistles') "
        "could be understood as an alteration in how these seminal reasons unfold after humanity's sin, "
        "thus describing a change in nature's operation."
    )

    print(f"Theologian: {theologian}\n")
    print("Explanation:")
    # Use textwrap to format the explanation neatly
    print(textwrap.fill(explanation, width=80))

if __name__ == "__main__":
    theological_answer()
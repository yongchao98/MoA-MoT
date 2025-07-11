import textwrap

def explain_biological_finding():
    """
    This function explains the answer to the biological question based on
    findings from a specific scientific paper.
    """
    # Key information from the user's question
    organism = "Microcystis aeruginosa CAAT 2005-3"
    stressor = "High temperature exposure (29ºC)"
    question = "Which antioxidants are initially activated to counteract oxidative stress?"

    # Findings from the relevant scientific paper (Lopes, C.S., et al., 2020)
    primary_finding = ("In response to heat stress, the cells selectively activated the "
                       "enzymatic antioxidants, specifically superoxide dismutase (SOD) "
                       "and catalase (CAT).")

    secondary_finding = ("In contrast, the study observed no initial changes in the "
                         "content of liposoluble antioxidants (like tocopherols) or "
                         "photosynthetic pigments.")

    conclusion = "Based on these findings, the correct answer is C: Enzymatic antioxidants."

    # Print the explanation in a clear format
    print("Step-by-Step Analysis of the Scientific Question:")
    print("-" * 50)
    print(f"Organism of Interest: {organism}")
    print(f"Environmental Stressor: {stressor}")
    print(f"Question: {question}")
    print("-" * 50)
    print("Conclusion from Published Research:")
    print("\n".join(textwrap.wrap(primary_finding, width=70)))
    print("\n")
    print("\n".join(textwrap.wrap(secondary_finding, width=70)))
    print("-" * 50)
    print(f"Final Answer: {conclusion}")


explain_biological_finding()
import textwrap

def analyze_aphid_metabolism():
    """
    Analyzes the provided statements about aphid metabolism based on the experimental description.
    """
    # Define the premises from the experiment description.
    premises = [
        "CA (watermelon-adapted) aphid does well on a high-raffinose diet.",
        "MA (cotton-adapted) aphid does well on a sucrose-only diet.",
        "Raffinose is a Raffinose Family Oligosaccharide (RFO), metabolized by the enzyme galactosidase.",
        "Host transfer moves CA to cotton (implied low-raffinose, high-sucrose) and MA to watermelon (implied high-raffinose)."
    ]

    print("### Analysis of the Statements ###\n")
    print("Based on the following premises:")
    for p in premises:
        print(f"- {p}")
    print("\n" + "="*40 + "\n")

    # Analyze each choice
    choices = {
        'A': "CA biotypes have an enhanced ability to metabolize RFOs than MA biotypes.",
        'B': "CA preferred raffinose-rich diet, whereas MA preferred sucrose-rich diet.",
        'C': "Upon the host transfer, the CA biotypes showed decreased galactosidase activity, owing to lower raffinose levels in cotton.",
        'D': "Upon the host transfer, the CA biotypes showed decreased galactosidase activity, owing to higher glucose levels in cotton.",
        'E': "Upon the host transfer, the MA biotypes showed increased galactosidase activity, owing to higher raffinose levels in the watermelon."
    }

    # Evaluation
    print("Choice A: " + choices['A'])
    analysis_a = "This is strongly supported. CA's success on a raffinose-rich diet implies it has the necessary metabolic machinery (like galactosidase), which MA, preferring sucrose, would likely have in lower amounts. Conclusion: LIKELY TRUE."
    print(textwrap.fill(f"   Analysis: {analysis_a}\n", width=80))

    print("Choice B: " + choices['B'])
    analysis_b = "The phrase 'did well' is a strong indicator of metabolic suitability, which can be interpreted as a preference in this biological context. Conclusion: LIKELY TRUE."
    print(textwrap.fill(f"   Analysis: {analysis_b}\n", width=80))

    print("Choice C: " + choices['C'])
    analysis_c = "CA aphids are moved to cotton, which is the host for the sucrose-specialist MA. This implies cotton has low raffinose levels. Galactosidase metabolizes raffinose, so its activity would decrease due to the lack of its specific substrate. Conclusion: LIKELY TRUE."
    print(textwrap.fill(f"   Analysis: {analysis_c}\n", width=80))
    
    print("Choice D: " + choices['D'])
    analysis_d = "This statement suggests that decreased galactosidase activity is caused by high glucose (from sucrose breakdown in cotton). While biologically plausible (a mechanism called catabolite repression), the text does not provide evidence for this specific mechanism. The more direct cause, as stated in C, is the lack of raffinose. Attributing the effect to glucose levels is an assumption not supported by the text. Therefore, this statement is not necessarily true. Conclusion: LIKELY NOT TRUE."
    print(textwrap.fill(f"   Analysis: {analysis_d}\n", width=80))

    print("Choice E: " + choices['E'])
    analysis_e = "MA aphids are moved to watermelon, a raffinose-rich environment. The presence of a new substrate (raffinose) would induce the production and activity of its corresponding enzyme (galactosidase) to metabolize it. Conclusion: LIKELY TRUE."
    print(textwrap.fill(f"   Analysis: {analysis_e}\n", width=80))
    
    print("="*40)
    print("Final Conclusion: Statement D makes an assumption about a specific biological mechanism (repression by glucose) that is not supported by the provided text. The other statements are direct and logical conclusions from the information given.")

# Run the analysis
analyze_aphid_metabolism()
<<<D>>>
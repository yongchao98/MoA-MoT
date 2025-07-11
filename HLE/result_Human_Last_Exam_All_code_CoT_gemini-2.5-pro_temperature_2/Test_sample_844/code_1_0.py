def analyze_aphid_statements():
    """
    Analyzes statements about aphid biotypes based on provided experimental data.
    """

    # --- Step 1: Establish Facts from the Text ---
    facts = {
        'CA_biotype': {
            'adaptation': 'watermelon',
            'diet_preference': 'raffinose-rich',
            'metabolizes': 'raffinose',
            'enzyme': 'galactosidase'
        },
        'MA_biotype': {
            'adaptation': 'cotton',
            'diet_preference': 'sucrose-rich',
            'metabolizes': 'sucrose'
        },
        'hosts': {
            'watermelon': 'raffinose-rich',
            'cotton': 'sucrose-rich'
        }
    }

    print("Analyzing each statement based on the experimental context...")
    print("-" * 60)

    # --- Step 2: Evaluate Each Statement Logically ---

    # A. CA biotypes have an enhanced ability to metabolize RFOs than MA biotypes.
    # Raffinose is an RFO. CA thrives on it, MA thrives on sucrose.
    analysis_A = (
        "TRUE. The text states that the CA biotype thrives on a raffinose-rich diet, while the MA biotype "
        "thrives on sucrose. This directly implies that CA has a superior ability to metabolize "
        "raffinose (an RFO) compared to MA."
    )
    print("Evaluation of A:", analysis_A)

    # B. CA preferred raffinose-rich diet, whereas MA preferred sucrose-rich diet.
    # This is a direct summary of the first sentence.
    analysis_B = (
        "TRUE. This statement is a direct summary of the information given in the first sentence of the "
        "prompt where it describes what diet each biotype 'did well on'."
    )
    print("Evaluation of B:", analysis_B)

    # C. Upon the host transfer, the CA biotypes showed decreased galactosidase activity, owing to lower raffinose levels in cotton.
    # CA moves from watermelon (high raffinose) to cotton (low raffinose).
    # Less substrate (raffinose) leads to lower activity of the corresponding enzyme (galactosidase).
    analysis_C = (
        "TRUE. The CA biotype is transferred from its native host watermelon (inferred to be raffinose-rich) "
        "to cotton (inferred to be low in raffinose). A decrease in the substrate (raffinose) "
        "would logically lead to decreased activity of its metabolizing enzyme, galactosidase."
    )
    print("Evaluation of C:", analysis_C)
    
    # D. Upon the host transfer, the CA biotypes showed decreased galactosidase activity, owing to higher glucose levels in cotton.
    # This suggests a mechanism called catabolite repression. While plausible, another statement is definitively false.
    analysis_D = (
        "PLAUSIBLE. This describes a valid biological mechanism (catabolite repression) that could cause "
        "decreased galactosidase activity. However, we are looking for the statement that is 'not true'. "
        "This statement describes a possible reality, unlike statement E, which contradicts the core premise."
    )
    print("Evaluation of D:", analysis_D)

    # E. Upon the host transfer, the MA biotypes showed increased galactosidase activity, owing to higher raffinose levels in the watermelon.
    # MA is adapted to sucrose, not raffinose. It's unlikely to possess the ability to *increase* enzyme activity for a sugar it is not adapted to.
    analysis_E = (
        "NOT TRUE. The MA biotype is described as adapted to sucrose and does well on a sucrose-only diet, "
        "implying it is poorly equipped to metabolize raffinose. It is highly unlikely that an organism "
        "would show *increased* activity of an enzyme for a compound it is not adapted to process. "
        "Instead, it would be expected to perform poorly."
    )
    print("Evaluation of E:", analysis_E)
    
    print("-" * 60)
    print("Conclusion: Statement E contradicts the established facts about the MA biotype's adaptation.")

# Run the analysis
analyze_aphid_statements()
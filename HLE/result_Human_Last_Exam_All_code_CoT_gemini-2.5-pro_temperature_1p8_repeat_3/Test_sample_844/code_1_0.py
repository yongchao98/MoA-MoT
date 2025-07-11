import textwrap

def analyze_aphid_statements():
    """
    Analyzes statements about aphid biotypes based on a given text.
    """
    # 1. Deconstruct the text into a knowledge base
    knowledge = {
        "CA_biotype": {
            "name": "CA",
            "adaptation": "watermelon",
            "diet": "raffinose-rich (sucrose:raffinose 3:8)",
            "sugar": "raffinose",
            "enzyme": "galactosidase"
        },
        "MA_biotype": {
            "name": "MA",
            "adaptation": "cotton",
            "diet": "sucrose-only",
            "sugar": "sucrose",
            "enzyme": "sucrase (implies low/no galactosidase capability)"
        },
        "inferences": {
            "watermelon": "raffinose-rich",
            "cotton": "raffinose-poor, sucrose-rich"
        }
    }

    print("Analyzing the statements based on the provided text...\n")

    # 2. Evaluate each statement
    analysis = []

    # Statement A: CA biotypes have an enhanced ability to metabolize RFOs than MA biotypes.
    # Raffinose is an RFO. CA thrives on raffinose, MA thrives on sucrose.
    # This statement is consistent with the text.
    conclusion_A = (
        "TRUE. The text states the CA biotype thrives on a raffinose-rich diet, while the MA "
        "biotype thrives on a sucrose-only diet. Raffinose is a Raffinose Family Oligosaccharide (RFO). "
        "Therefore, the CA biotype is better adapted to metabolize RFOs."
    )
    analysis.append(('A', True, conclusion_A))


    # Statement B: CA preferred raffinose-rich diet, whereas MA preferred sucrose-rich diet.
    # "Did well on" implies preference or adaptation.
    # This statement is consistent with the text.
    conclusion_B = (
        "TRUE. This is a direct summary of the first sentence of the text, where the performance "
        "of each biotype on different artificial diets is described."
    )
    analysis.append(('B', True, conclusion_B))


    # Statement C: Upon the host transfer, the CA biotypes showed decreased galactosidase activity, owing to lower raffinose levels in cotton.
    # CA moves to cotton, which is inferred to be raffinose-poor. Less substrate (raffinose) means less need for the enzyme (galactosidase).
    # This is a sound biological conclusion.
    conclusion_C = (
        "TRUE. The CA biotype is adapted to raffinose (likely high in watermelon). When transferred to cotton, "
        "which is sucrose-rich and thus raffinose-poor, the substrate for galactosidase is scarce. "
        "In biological systems, enzyme activity often decreases in the absence of its substrate."
    )
    analysis.append(('C', True, conclusion_C))


    # Statement D: Upon the host transfer, the CA biotypes showed decreased galactosidase activity, owing to higher glucose levels in cotton.
    # This is similar to C but offers a different, though plausible, reason. High glucose (from sucrose) can inhibit enzymes for other sugars (catabolite repression).
    # This is a plausible biological mechanism. It's hard to deem it "not true".
    conclusion_D = (
        "LIKELY TRUE. This describes the same outcome as C but gives a different valid biological reason. "
        "High levels of a preferred sugar like glucose (from sucrose in cotton) can cause catabolite repression, "
        "inhibiting the production of enzymes for less-preferred sugars. This is a plausible explanation."
    )
    analysis.append(('D', True, conclusion_D))


    # Statement E: Upon the host transfer, the MA biotypes showed increased galactosidase activity, owing to higher raffinose levels in the watermelon.
    # The MA biotype is a sucrose specialist. It is not adapted to metabolize raffinose.
    # It is highly unlikely that an unadapted, specialized organism would simply increase enzyme activity for a novel food source.
    # This statement contradicts the established specialization of the MA biotype.
    conclusion_E = (
        "NOT TRUE. The text establishes that the MA biotype is a cotton/sucrose specialist. "
        "It is not adapted to metabolize raffinose and likely has very low or inactive genes for galactosidase. "
        "Transferring it to a raffinose-rich environment would likely result in poor survival, not an *increase* "
        "in the activity of an enzyme it is not specialized to produce. This contradicts its fundamental adaptation."
    )
    analysis.append(('E', False, conclusion_E))


    # 3. Print the analysis and find the "not true" statement
    not_true_statement = None
    for choice, is_true, reason in analysis:
        print(f"Choice {choice}:")
        # textwrap helps format long strings neatly
        print(textwrap.fill(reason, width=80))
        print("-" * 20)
        if not is_true:
            not_true_statement = choice

    print(f"\nThe statement that is not true is: {not_true_statement}")
    print(f"\n<<<E>>>")


if __name__ == '__main__':
    analyze_aphid_statements()
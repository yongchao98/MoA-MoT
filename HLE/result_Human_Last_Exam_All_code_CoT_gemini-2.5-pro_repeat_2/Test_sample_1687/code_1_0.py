import operator

def diagnose_case():
    """
    This function analyzes clinical findings from a case vignette to determine the most likely diagnosis.
    It uses a scoring system where clinical findings are weighted based on their relevance to each diagnosis.
    """

    # Step 1: Define the key clinical findings present in the patient's case.
    # The absence of a finding can also be a key piece of information.
    case_findings = {
        "difficult colonoscopy": True,
        "LUQ pain": True,
        "left shoulder pain (Kehr's sign)": True,
        "acute massive blood loss (dropping Hgb)": True,
        "hemodynamic instability (shock)": True,
        "peritoneal signs": True,
        "no polypectomy performed": True  # A crucial negative finding
    }

    # Step 2: Define the criteria and their associated scores for each diagnosis.
    # High scores for strongly suggestive findings, negative scores for contradictory findings.
    diagnosis_criteria = {
        "A. Colonic perforation": {
            "peritoneal signs": 2,
            "LUQ pain": 1, # Pain is present, but location is less typical for colon perforation.
            "hemodynamic instability (shock)": 1,
            "left shoulder pain (Kehr's sign)": 0 # Not a typical sign.
        },
        "B. Lower GI bleeding": {
            "acute massive blood loss (dropping Hgb)": 2,
            "hemodynamic instability (shock)": 1,
            "LUQ pain": -1, # Pain is typically lower abdominal, not LUQ.
            "left shoulder pain (Kehr's sign)": -2 # Not a sign of lower GI bleeding.
        },
        "C. Splenic laceration": {
            "LUQ pain": 2, # Classic location for spleen injury.
            "left shoulder pain (Kehr's sign)": 3, # Highly specific sign for splenic injury.
            "acute massive blood loss (dropping Hgb)": 2, # Hallmark of splenic rupture.
            "hemodynamic instability (shock)": 2, # Expected with massive hemorrhage.
            "difficult colonoscopy": 1 # Known risk factor (traction on splenocolic ligament).
        },
        "D. Postpolypectomy syndrome": {
            "no polypectomy performed": -10, # This finding definitively rules out the diagnosis.
            "peritoneal signs": 1,
            "LUQ pain": 1
        }
    }

    # Step 3: Calculate the score for each diagnosis and print the breakdown.
    scores = {}
    print("Calculating likelihood scores for each diagnosis:\n")

    for diagnosis, criteria in diagnosis_criteria.items():
        total_score = 0
        equation_parts = []
        for finding, score in criteria.items():
            # Check if the finding is present in the case
            if case_findings.get(finding, False):
                total_score += score
                equation_parts.append(f"{score} (for {finding})")
        
        scores[diagnosis] = total_score
        
        # Print the equation for each diagnosis
        equation_str = " + ".join(equation_parts)
        if not equation_str:
            equation_str = "0"
            
        print(f"Score for {diagnosis}:\n  {equation_str} = {total_score}\n")

    # Step 4: Determine the most likely diagnosis
    most_likely_diagnosis = max(scores.items(), key=operator.itemgetter(1))[0]

    print("--------------------------------------------------")
    print(f"The most likely diagnosis is the one with the highest score.")
    print(f"Conclusion: {most_likely_diagnosis}")
    print("--------------------------------------------------")


if __name__ == "__main__":
    diagnose_case()
<<<C>>>
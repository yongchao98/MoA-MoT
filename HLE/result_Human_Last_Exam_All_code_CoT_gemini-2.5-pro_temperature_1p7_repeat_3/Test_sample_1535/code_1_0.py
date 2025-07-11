import textwrap

def solve_medical_case():
    """
    This function explains the reasoning behind the diagnosis based on the provided clinical vignette.
    """
    
    # Patient Presentation
    patient_findings = {
        "History": "muscle weakness, myalgia, fatigue",
        "Physical Exam": "periorbital recession and erythema"
    }
    
    print("Step 1: Identify key signs and symptoms from the patient's presentation.")
    print(f"- Muscle-related symptoms: {patient_findings['History']}")
    print(f"- Key skin finding on exam: {patient_findings['Physical Exam']}\n")
    
    print("Step 2: Connect the key findings to a likely diagnosis.")
    explanation_diagnosis = ("The combination of muscle inflammation symptoms (myalgia, weakness) and a "
                             "characteristic skin rash is highly suggestive of Dermatomyositis.")
    print(textwrap.fill(explanation_diagnosis, width=80))
    print("")

    print("Step 3: Define the specific skin finding and its location.")
    explanation_rash = ("The term 'periorbital erythema' describes redness (erythema) around the eyes (periorbital). "
                        "This finding in Dermatomyositis is known as a 'heliotrope rash'. "
                        "This rash is pathognomonic (specifically characteristic) of the disease and is located on the eyelids.")
    print(textwrap.fill(explanation_rash, width=80))
    print("")
    
    print("Step 4: Conclude based on the evidence.")
    conclusion = ("Given that the patient has 'periorbital erythema' mentioned in the physical exam, the anatomical region "
                  "expected to have this rash is the eyelids.")
    print(textwrap.fill(conclusion, width=80))
    print("-" * 30)

solve_medical_case()
print("Final Answer based on the clinical description:")
print("C")
<<<C>>>
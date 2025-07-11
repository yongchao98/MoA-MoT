def diagnose_newborn():
    """
    Analyzes clinical findings to determine the most likely anatomical defect.
    """

    # Clinical findings from the prompt
    patient_findings = {
        "weight_lb": 12,
        "weight_oz": 1,
        "o2_sat_percent": 89,
        "radiology": "fluid-filled density in the left lung",
        "physical_exam": "micrognathia"
    }

    # Answer choices
    answer_choices = {
        "A": "Pulmonary hypoplasia",
        "B": "Pleuroperitoneal membrane defect",
        "C": "Ventral foregut budding defect",
        "D": "Maternal diabetes",
        "E": "Fistula"
    }

    # Explanation based on clinical reasoning
    analysis = {
        "A": "Causes respiratory distress but is often a secondary finding to lung compression.",
        "B": "Explains all findings: defect leads to diaphragmatic hernia (left lung density), which causes lung compression (respiratory distress/hypoplasia). Associated with micrognathia and maternal diabetes risk factors.",
        "C": "Leads to tracheoesophageal fistula; doesn't explain the specific lung density or micrognathia as well.",
        "D": "This is a maternal risk factor for the baby's condition (macrosomia, CDH), not the baby's anatomical defect.",
        "E": "Too general. A pleuroperitoneal membrane defect is a specific pathology that fits perfectly."
    }

    # Based on the synthesis, B is the most comprehensive explanation.
    most_likely_defect_key = "B"
    most_likely_defect_desc = answer_choices[most_likely_defect_key]

    print("Patient's Clinical Data:")
    print(f"- Weight: {patient_findings['weight_lb']} lb {patient_findings['weight_oz']} oz")
    print(f"- Oxygen Saturation: {patient_findings['o2_sat_percent']}%")
    print(f"- Key Findings: Left lung density and micrognathia.")
    print("\n--- Analysis ---")
    print("Evaluating the connection between findings and potential defects:")
    print(f"The combination of a left-sided chest mass, respiratory distress, and micrognathia is characteristic of a Congenital Diaphragmatic Hernia (CDH).")
    print("The most common cause of CDH is a failure of the pleuroperitoneal membrane to close during development.")
    print("\n--- Conclusion ---")
    print(f"The most likely anatomical defect is: ({most_likely_defect_key}) {most_likely_defect_desc}")

diagnose_newborn()
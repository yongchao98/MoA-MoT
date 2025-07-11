def analyze_medical_case():
    """
    This script analyzes the clinical findings from the case study to determine the most likely diagnosis.
    It evaluates each possibility based on the key symptoms and lab results.
    """

    # Key findings extracted from the case report
    initial_hemoglobin = 11.7
    post_bleed_hemoglobin = 6.5
    pain_location = "Left Upper Quadrant (LUQ) and Left Shoulder"
    procedure_details = "Difficult colonoscopy, but NO polypectomy performed"
    signs_of_internal_bleeding = ["Worsening abdominal distension", "Peritoneal signs", "Tachycardia and hypotension (shock)"]

    print("Analyzing the patient's clinical presentation step-by-step:\n")

    # Analysis of Diagnosis A: Colonic Perforation
    print("--- Evaluating Diagnosis A: Colonic Perforation ---")
    print("This is less likely because the massive blood loss (hemoglobin drop from {} to {}) and specific pain in the LUQ and left shoulder are not the most typical signs. The presentation points more strongly to a solid organ injury.\n".format(initial_hemoglobin, post_bleed_hemoglobin))

    # Analysis of Diagnosis B: Lower GI Bleeding
    print("--- Evaluating Diagnosis B: Lower GI Bleeding ---")
    print("This is unlikely. The patient's signs of abdominal distension and peritoneal irritation indicate bleeding INTO the abdominal cavity (hemoperitoneum), not bleeding OUT through the rectum.\n")

    # Analysis of Diagnosis C: Splenic Laceration
    print("--- Evaluating Diagnosis C: Splenic Laceration ---")
    print("This is the most likely diagnosis. It perfectly aligns with all the key findings:")
    print(f"*   Location: The spleen is in the LUQ, matching the '{pain_location}'.")
    print("*   Kehr's Sign: Left shoulder pain is a classic sign of splenic injury causing diaphragmatic irritation.")
    print(f"*   Massive Hemorrhage: The spleen is highly vascular, explaining the severe hemoglobin drop from {initial_hemoglobin} g/dL to {post_bleed_hemoglobin} g/dL and the resulting shock.")
    print(f"*   Hemoperitoneum: The signs of {', '.join(signs_of_internal_bleeding)} are explained by blood filling the abdomen.\n")

    # Analysis of Diagnosis D: Postpolypectomy Syndrome
    print("--- Evaluating Diagnosis D: Postpolypectomy Syndrome ---")
    print(f"This is ruled out because the case explicitly states: '{procedure_details}'.\n")

    # Final Conclusion
    print("Conclusion: The combination of LUQ pain, referred left shoulder pain, and profound hypovolemic shock after a difficult colonoscopy makes Splenic Laceration the most probable diagnosis.")


analyze_medical_case()
<<<C>>>
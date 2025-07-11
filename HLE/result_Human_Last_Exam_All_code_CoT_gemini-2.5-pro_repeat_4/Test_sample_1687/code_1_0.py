import textwrap

def solve_medical_case():
    """
    Analyzes the clinical case and determines the most likely diagnosis.
    """

    # Key numerical values extracted from the case description
    patient_age = 65
    initial_hemoglobin = 11.7
    post_deterioration_hemoglobin = 6.5
    final_hemoglobin = 7.1
    heart_rate = 105
    blood_pressure_systolic = 100
    blood_pressure_diastolic = 65

    # Step-by-step reasoning
    reasoning = f"""
    The diagnosis is determined by analyzing the unique combination of symptoms and lab results following a specific medical procedure.

    1.  **Procedure and Location of Pain:** The patient had a difficult colonoscopy and subsequently developed pain in the Left Upper Quadrant (LUQ) and referred pain to the left shoulder. The spleen is located in the LUQ. Pain in the left shoulder (Kehr's sign) is a classic sign of diaphragmatic irritation from splenic bleeding.

    2.  **Evidence of Severe Hemorrhage:** The patient's clinical course demonstrates a massive internal bleed.
        *   Her hemoglobin level dropped precipitously from {initial_hemoglobin} g/dL to {post_deterioration_hemoglobin} g/dL.
        *   She developed signs of hemorrhagic shock: tachycardia (heart rate {heart_rate} beats/min) and hypotension (blood pressure {blood_pressure_systolic}/{blood_pressure_diastolic} mm Hg).
        *   Physical signs included pallor and increasing abdominal distension with peritoneal signs, consistent with a large volume of blood in the abdomen (hemoperitoneum).

    3.  **Evaluating the Options:**
        *   (A) Colonic perforation is less likely to cause this degree of rapid blood loss and does not typically present with referred left shoulder pain.
        *   (B) Lower GI bleeding would manifest as blood in the stool, which was not reported. The symptoms point to an intra-abdominal, not intraluminal, bleed.
        *   (D) Postpolypectomy syndrome is ruled out as no polypectomy was performed.
        *   (C) Splenic laceration is a rare complication of colonoscopy that perfectly explains the entire clinical picture: the specific location of pain, the classic referred pain pattern, and the profound hemorrhagic shock.

    Conclusion: The constellation of findings strongly points to a splenic injury.
    """

    print(textwrap.dedent(reasoning).strip())

    # Final Answer
    final_answer = "C"
    print(f"\nFinal Answer:\n<<<{final_answer}>>>")

solve_medical_case()
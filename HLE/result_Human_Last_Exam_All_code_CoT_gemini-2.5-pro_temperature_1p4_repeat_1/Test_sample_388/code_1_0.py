def generate_counseling_recommendation():
    """
    Analyzes patient medications to provide a key counseling point.
    """
    # Patient's prescriptions and OTC medication
    otc_med = "Excedrin (contains Aspirin)"
    prescription_1 = "Atorvastatin 20mg"
    prescription_2 = "Junel Fe 1.5/30mg"
    prescription_3 = "Fluoxetine 20mg"

    # The core recommendation based on the interaction between Aspirin and Fluoxetine
    recommendation = f"""
Pharmacist Counseling Recommendation:

The most important recommendation is based on a potential interaction between the Excedrin you took and one of your prescriptions.

1.  **The Interaction:** Your prescription for {prescription_3} can have a serious interaction with the aspirin component found in Excedrin.

2.  **The Risk:** Taking these two medications together significantly increases your risk of bleeding, especially stomach bleeding.

3.  **The Recommendation:** For future headaches, you should avoid products containing aspirin or other NSAIDs (like ibuprofen). A safer choice for you would be a product that contains only acetaminophen, such as Tylenol. It is effective for pain relief without increasing this specific bleeding risk.

While your other prescriptions, {prescription_1} and {prescription_2}, do not have this specific interaction, it is always good practice to discuss any new over-the-counter medication use with your pharmacist or doctor.
"""
    print(recommendation)

generate_counseling_recommendation()
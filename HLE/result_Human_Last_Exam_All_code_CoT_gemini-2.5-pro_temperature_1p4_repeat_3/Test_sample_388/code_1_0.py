def pharmacist_counseling():
    """
    This function formulates and prints a pharmacist's counseling recommendation
    based on a potential drug interaction.
    """
    prescription_1 = "Atorvastatin 20mg"
    prescription_2 = "Junel Fe 1.5/30mg"
    prescription_3 = "Fluoxetine 20mg"
    otc_drug = "Excedrin"
    interacting_ingredient_in_otc = "Aspirin (an NSAID)"
    interacting_prescription = "Fluoxetine (an SSRI)"

    recommendation = f"""
Based on the medicines you are taking, there is an important counseling point regarding a potential drug interaction.

Interaction Identified:
- Prescription: {interacting_prescription}
- Over-the-counter medicine: {otc_drug}, which contains {interacting_ingredient_in_otc}

Counseling Recommendation:
"I see you are picking up your prescription for {prescription_3.split(' ')[0]} and that you recently used {otc_drug} for a headache. It's important to be aware that taking these two together can increase the risk of stomach bleeding.

For future headaches, a safer alternative would be a product containing only acetaminophen, like Tylenol, as it does not have this same interaction risk. Please watch for any signs of bleeding, such as unusual bruising, dark or tarry stools, or stomach pain, and let your doctor know if they occur."
"""
    print(recommendation)

pharmacist_counseling()
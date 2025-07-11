def provide_counseling_recommendation():
    """
    Analyzes Allison's medications to provide a key counseling point.
    """
    
    # Medication details from the user's prompt
    prescription_1 = "Atorvastatin"
    dosage_1 = 20  # mg

    prescription_2 = "Junel Fe"
    dosage_2_part_1 = 1.5  # mg
    dosage_2_part_2 = 30  # mg, but this is usually mcg, so let's stick to the numbers given

    prescription_3 = "Fluoxetine"
    dosage_3 = 20  # mg

    otc_med = "Excedrin (which contains Aspirin, an NSAID)"

    # The final recommendation text
    recommendation = f"""
Based on Allison's medications, the most critical counseling recommendation involves a drug interaction.

**Identified Interaction:**
The combination of {prescription_3} (an SSRI antidepressant) and the aspirin in {otc_med} (an NSAID) significantly increases the risk of bleeding, especially stomach bleeding.

**Pharmacist's Recommendation:**
Allison should be advised to AVOID taking Excedrin or other products containing aspirin or NSAIDs (like ibuprofen) while she is on {prescription_3}.

**Safer Alternative for Headache:**
A safer option would be a product containing only acetaminophen (e.g., Tylenol).

**Summary of Prescription Dosages:**
- {prescription_1}: {dosage_1}mg
- {prescription_2}: {dosage_2_part_1}/{dosage_2_part_2}mg
- {prescription_3}: {dosage_3}mg
"""
    
    print(recommendation)

provide_counseling_recommendation()
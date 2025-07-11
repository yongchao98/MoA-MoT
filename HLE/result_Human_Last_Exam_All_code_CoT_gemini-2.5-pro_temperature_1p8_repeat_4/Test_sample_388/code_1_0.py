def provide_pharmacist_counseling():
    """
    Analyzes a patient's medication profile to generate key counseling recommendations.
    """
    # Define the patient's medications and their relevant doses/classes.
    fluoxetine = {"name": "Fluoxetine", "dose": "20mg", "class": "SSRI"}
    junel_fe = {"name": "Junel Fe 1.5/30mg", "class": "Oral Contraceptive"}
    atorvastatin = {"name": "Atorvastatin", "dose": "20mg"}
    otc_med = {"name": "Excedrin", "key_ingredient": "Aspirin"}

    print("Based on the information provided, here are the most important counseling recommendations:\n")

    # Point 1: The most critical interaction between Fluoxetine and Excedrin.
    print("1. Bleeding Risk with Fluoxetine and Excedrin:")
    print(f"   - Your new prescription for {fluoxetine['name']} ({fluoxetine['dose']}) can have a serious interaction with the {otc_med['key_ingredient']} found in the {otc_med['name']} you took for your headache.")
    print(f"   - Taking an SSRI like {fluoxetine['name']} with an NSAID like {otc_med['key_ingredient']} significantly increases the risk of bleeding, especially in the stomach.")
    print("   - RECOMMENDATION: For future headaches, you should use a product that only contains acetaminophen (like Tylenol) and avoid products with aspirin or ibuprofen.\n")

    # Point 2: The potential seriousness of a headache while on oral contraceptives.
    print("2. Headache Monitoring with Junel Fe:")
    print(f"   - It is also very important to be aware that for your prescription {junel_fe['name']}, a new or worsening headache can be a rare but serious side effect.")
    print("   - This could be a warning sign of a blood clot.")
    print("   - RECOMMENDATION: Since you had a headache today, please monitor it carefully. If it returns, is worse than usual, or is accompanied by vision changes or numbness, you must contact your doctor right away.\n")

    # Point 3: Standard counseling for Atorvastatin.
    print("3. General Counseling for Atorvastatin:")
    print(f"   - For your {atorvastatin['name']} ({atorvastatin['dose']}), a common counseling point is to be aware of muscle pain. If you notice any new or unexplained muscle pain or weakness, let your doctor know.")


# Execute the function to print the counseling points.
provide_pharmacist_counseling()
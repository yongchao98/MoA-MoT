def provide_pharmacist_counseling():
    """
    This script outlines a pharmacist's counseling points for Allison.
    """
    
    # Prescription information from the scenario
    medications = {
        "Atorvastatin": "20mg",
        "Junel Fe": "1.5/30mg",
        "Fluoxetine": "20mg"
    }
    
    print("Patient's Prescriptions:")
    for med, dosage in medications.items():
        print(f"- {med}: {dosage}")
        
    print("\n" + "="*50 + "\n")
    
    print("Pharmacist Counseling Recommendation:")
    
    recommendation = """
Hello Allison, I've got your prescriptions ready for you.

I overheard you mentioning that you had a headache. I see that you are taking Junel Fe, and it's important for you to know that new or worsening headaches can sometimes be a side effect of combination birth control pills like this one.

Because of this, I have an important question for you: Before your headache starts, do you ever experience any unusual visual symptoms, like seeing flashing lights, zig-zag lines, or blind spots? Or any numbness or tingling in your hand or face? This is known as an 'aura.'

This is a crucial safety question because having migraines with aura while taking a birth control pill with estrogen increases the risk of a stroke.

Regardless of whether you have these specific symptoms, since these headaches are a new problem for you, I strongly recommend that you contact the doctor who prescribed your Junel Fe to let them know. They need to be aware of this new symptom to ensure this is still the safest and best option for you. They might want to discuss other birth control options if these headaches are related to your medication.
"""
    
    print(recommendation)

if __name__ == "__main__":
    provide_pharmacist_counseling()
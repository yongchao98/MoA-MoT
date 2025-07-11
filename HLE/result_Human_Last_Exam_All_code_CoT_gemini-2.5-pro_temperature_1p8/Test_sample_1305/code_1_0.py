def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis by printing a step-by-step evaluation.
    """
    print("Step 1: Analyzing the patient's presentation.")
    print(" - History: 31-year-old woman with progressive shortness of breath, fatigue, and a history of faint childhood cyanosis. This suggests a congenital heart defect with long-term consequences.")
    print(" - Physical Exam: A systolic ejection murmur at the left upper sternal border (LUSB) which increases with inspiration. This is a classic finding for a right-sided heart lesion causing increased flow or obstruction, such as a pulmonic stenosis murmur or an atrial septal defect (ASD) flow murmur.")
    print(" - ECG: Right ventricular hypertrophy (RVH) and Left Axis Deviation (LAD). RVH is expected with a right-sided lesion. However, the LAD is a very specific clue, as most right-sided lesions cause Right Axis Deviation.")

    print("\nStep 2: Evaluating the most likely causes based on the unique findings.")
    print(" - The task is to find a diagnosis that explains a right-sided volume/pressure overload (murmur at LUSB increasing with inspiration, RVH) AND a left axis deviation.")
    print(" - Most congenital defects causing these right-sided findings (like a standard secundum ASD or pulmonic stenosis) are associated with Right Axis Deviation.")
    print(" - We must look for an exception.")

    print("\nStep 3: Pinpointing the correct diagnosis.")
    print(" - Atrial Septal Defect (ASD) is a strong candidate for the murmur type. An ASD causes a left-to-right shunt, increasing blood volume across the pulmonic valve, creating a flow murmur that matches the description.")
    print(" - While the most common type of ASD (ostium secundum) causes right axis deviation, a less common type called an 'ostium primum' ASD is classically associated with Left Axis Deviation.")
    print(" - An ostium primum ASD perfectly reconciles the findings: It causes right ventricular volume overload (leading to the described flow murmur and RVH) and simultaneously causes Left Axis Deviation due to the abnormal position of the defect affecting the heart's conduction system.")
    print(" - Therefore, an ASD (specifically, an ostium primum type) is the most likely diagnosis.")

    print("\nStep 4: Final Conclusion.")
    print(" - The most fitting answer choice that encompasses this diagnosis is D. Atrial septal defect.")

solve_medical_case()
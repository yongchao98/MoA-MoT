def diagnose_skin_condition():
    """
    This script analyzes clinical findings to determine the most likely diagnosis.
    It assigns scores to key features to quantify the evidence for Hidradenitis Suppurativa (HS).
    """

    # --- Patient's Key Clinical Findings ---
    # We assign a score based on the diagnostic importance of each finding for HS.
    findings = {
        "Presence of Purulent Nodules": 3,
        "Lesions in Multiple Intertriginous Areas (axillae, groin, inframammary)": 3,
        "Risk Factor: Smoking": 2,
        "Risk Factor: Obesity (BMI 39)": 2
    }

    print("Analyzing clinical evidence for Hidradenitis Suppurativa...")
    print("---------------------------------------------------------")

    total_score = 0
    equation_parts = []

    for finding, score in findings.items():
        print(f"Evidence found: '{finding}'. Diagnostic score contribution: {score}")
        total_score += score
        equation_parts.append(str(score))

    # --- Final Calculation and Conclusion ---
    print("\nCalculating the final diagnostic likelihood score...")
    
    # Constructing and printing the final equation as requested
    final_equation = " + ".join(equation_parts)
    print(f"Final Score Equation: {final_equation} = {total_score}")

    print("\nConclusion:")
    print("The high score, based on the presence of purulent nodules in characteristic locations")
    print("combined with major risk factors (smoking and obesity), strongly points to a")
    print("diagnosis of Hidradenitis Suppurativa.")


if __name__ == "__main__":
    diagnose_skin_condition()

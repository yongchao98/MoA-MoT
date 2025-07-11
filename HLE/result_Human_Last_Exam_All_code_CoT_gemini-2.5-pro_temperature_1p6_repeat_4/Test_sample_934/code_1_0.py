def analyze_esophageal_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by assigning illustrative scores to key findings.
    """
    print("Analyzing patient factors to determine the most likely diagnosis...")

    # Key patient findings
    heavy_smoking = {"name": "Heavy Smoking History (>20 years)", "points": 4}
    alcohol_use_disorder = {"name": "Alcohol Use Disorder", "points": 3}
    imaging_shows_thickening = {"name": "Imaging shows wall thickening/narrowing", "points": 2}
    normal_endoscopy = {"name": "Normal endoscopy with abnormal imaging (suggests submucosal mass)", "points": 5}

    findings = [heavy_smoking, alcohol_use_disorder, imaging_shows_thickening, normal_endoscopy]

    total_score = 0
    equation_numbers = []

    print("\nCalculating a likelihood score for Esophageal Squamous Cell Carcinoma (SCC):")
    for finding in findings:
        score = finding["points"]
        total_score += score
        equation_numbers.append(str(score))
        print(f"- Finding: '{finding['name']}'. Score: +{score}")

    # This fulfills the request to show each number in a final equation.
    final_equation_str = " + ".join(equation_numbers)
    print("\nIllustrative final score equation:")
    print(f"{final_equation_str} = {total_score}")

    print("\nConclusion:")
    print("The combination of major risk factors (smoking, alcohol) with the key finding of")
    print("wall thickening on imaging despite a normal endoscopy is highly suggestive of a")
    print("submucosal process. Among the choices, this presentation is most characteristic")
    print("of Esophageal Squamous Cell Carcinoma.")

# Run the analysis
analyze_esophageal_case()
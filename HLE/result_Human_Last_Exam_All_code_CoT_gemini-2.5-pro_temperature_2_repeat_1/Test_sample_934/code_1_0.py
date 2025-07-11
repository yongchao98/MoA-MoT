import sys

def calculate_risk_score():
    """
    Analyzes patient data to calculate a risk score for the most
    likely esophageal pathologies based on the provided vignette.
    """

    patient_findings = {
        "Heavy Smoking (>40 pack-years)": {
            "scc_points": 3,
            "adenocarcinoma_points": 1
        },
        "Alcohol Use Disorder": {
            "scc_points": 3,
            "adenocarcinoma_points": 0
        },
        "Odynophagia (Pain with swallowing)": {
            "scc_points": 2,
            "adenocarcinoma_points": 1
        },
        "Imaging (Wall thickening/narrowing)": {
            "scc_points": 2,
            "adenocarcinoma_points": 2
        },
        "Lack of GERD history (implied)": {
            "scc_points": 0,
            "adenocarcinoma_points": -2 # Strong negative factor for adenocarcinoma
        }
    }

    scc_score = 0
    adeno_score = 0
    
    print("Risk Score Calculation for Esophageal Squamous Cell Carcinoma (SCC):")
    scc_equation_parts = []
    for finding, points in patient_findings.items():
        point_value = points.get("scc_points", 0)
        if point_value != 0:
            scc_score += point_value
            # Print each number for the equation as requested
            print(f"Finding: '{finding}', Score: {point_value}")
            scc_equation_parts.append(str(point_value))
    print(f"Final SCC Score = {' + '.join(scc_equation_parts)} = {scc_score}\n")


    print("Risk Score Calculation for Esophageal Adenocarcinoma:")
    adeno_equation_parts = []
    for finding, points in patient_findings.items():
        point_value = points.get("adenocarcinoma_points", 0)
        if point_value != 0:
            adeno_score += point_value
            # Print each number for the equation as requested
            print(f"Finding: '{finding}', Score: {point_value}")
            # Format for addition
            if point_value > 0:
              adeno_equation_parts.append(str(point_value))
            else:
              adeno_equation_parts.append(f"({point_value})")

    print(f"Final Adenocarcinoma Score = {' + '.join(adeno_equation_parts)} = {adeno_score}\n")

    if scc_score > adeno_score:
        print("Conclusion: The risk profile strongly favors Esophageal Squamous Cell Carcinoma.")
    else:
        print("Conclusion: The risk profile favors Esophageal Adenocarcinoma.")
        
    # The other options are significantly less likely.
    # Herpes/Strep esophagitis would show ulcers/plaques on endoscopy.
    # GERD does not typically cause significant wall thickening with a normal-appearing mucosa.

calculate_risk_score()
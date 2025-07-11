import sys
# This script is for educational and illustrative purposes only and does not constitute medical advice.
# A real medical decision should be made by a qualified healthcare professional.

def recommend_htn_medications():
    """
    Recommends a 3-drug hypertension regimen based on patient profile and contraindications.
    """
    # Medications the patient (MM) is unable to take
    excluded_medications = {
        'eplerenone', 'finerenone', 'hydrochlorothiazide', 'indapamide',
        'bumetanide', 'furosemide', 'torsemide', 'metolazone', 'verapamil'
    }

    # Guideline-based medication options for resistant hypertension, with preferred agents listed first.
    medication_options = {
        "Calcium Channel Blocker": ['Amlodipine', 'Nifedipine', 'Felodipine'],
        "RAAS Inhibitor (ACEi)": ['Lisinopril', 'Ramipril', 'Benazepril'],
        "Thiazide-like Diuretic": ['Chlorthalidone', 'Hydrochlorothiazide', 'Indapamide']
    }

    recommendations = []
    
    # Logic to select one drug from each core class for a 3-drug regimen
    
    # 1. Select a Calcium Channel Blocker (first-line for African American patients)
    for drug in medication_options["Calcium Channel Blocker"]:
        if drug.lower() not in excluded_medications:
            recommendations.append(drug)
            break
            
    # 2. Select a RAAS Inhibitor
    for drug in medication_options["RAAS Inhibitor (ACEi)"]:
        if drug.lower() not in excluded_medications:
            recommendations.append(drug)
            break
            
    # 3. Select a Thiazide-like Diuretic (cornerstone of HTN therapy)
    for drug in medication_options["Thiazide-like Diuretic"]:
        if drug.lower() not in excluded_medications:
            recommendations.append(drug)
            break

    # Print the final recommendations
    if len(recommendations) == 3:
        print("Based on the patient's profile and restrictions, the following 3 medications are recommended to maximize hypertension treatment:")
        print(f"1. {recommendations[0]}")
        print(f"2. {recommendations[1]}")
        print(f"3. {recommendations[2]}")
        # The line below is necessary for the final answer extraction.
        sys.stdout.write(f'<<<{recommendations[0]}, {recommendations[1]}, {recommendations[2]}>>>')
    else:
        print("Could not generate a full 3-drug regimen from the available options.")

if __name__ == '__main__':
    recommend_htn_medications()
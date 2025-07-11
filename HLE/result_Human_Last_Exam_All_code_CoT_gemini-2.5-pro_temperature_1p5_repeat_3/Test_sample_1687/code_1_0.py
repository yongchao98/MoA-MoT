import sys

def solve_diagnosis():
    """
    Analyzes clinical findings using a scoring system to determine the most likely diagnosis.
    """
    # A: Colonic perforation, B: Lower GI bleeding, C: Splenic laceration, D: Postpolypectomy syndrome
    scores = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
    equations = {
        'A': '0', 'B': '0', 'C': '0', 'D': '0'
    }
    
    print("Analyzing clinical findings to determine the most likely diagnosis...")
    print("-" * 50)

    # --- Finding 1: Difficult colonoscopy ---
    # This increases the risk for iatrogenic injury.
    scores['A'] += 1
    scores['C'] += 1
    equations['A'] += ' + 1'
    equations['C'] += ' + 1'
    print("Finding 1: 'Slightly difficult' colonoscopy")
    print("Reasoning: A difficult procedure increases the risk of trauma to the colon (perforation) or adjacent structures (spleen).")
    print("--> Score Update: +1 for Colonic perforation, +1 for Splenic laceration.\n")
    
    # --- Finding 2: No polypectomy performed ---
    # This is a critical piece of negative information.
    scores['D'] -= 100 # This effectively rules out Postpolypectomy syndrome.
    equations['D'] += ' - 100'
    print("Finding 2: 'No polypectomy was performed'")
    print("Reasoning: This finding rules out postpolypectomy syndrome as a diagnosis.")
    print("--> Score Update: -100 for Postpolypectomy syndrome.\n")

    # --- Finding 3: Left upper quadrant (LUQ) and left-shoulder pain ---
    # This combination is highly specific.
    scores['A'] += 1 # A perforation at the splenic flexure could cause LUQ pain.
    scores['C'] += 5 # This is the classic location for splenic pain and referred pain (Kehr's sign).
    equations['A'] += ' + 1'
    equations['C'] += ' + 5'
    print("Finding 3: LUQ pain combined with left-shoulder discomfort")
    print("Reasoning: This is a classic presentation (Kehr's sign) for splenic injury, where blood irritates the diaphragm.")
    print("--> Score Update: +5 for Splenic laceration, +1 for Colonic perforation.\n")

    # --- Finding 4: Rapid drop in hemoglobin and hemodynamic instability ---
    # This indicates severe, active bleeding.
    scores['A'] += 2 # Perforation can cause bleeding, but a drop this fast suggests major vessel involvement.
    scores['B'] -= 2 # Intraluminal GI bleeding would likely present with blood per rectum, not just pain and shock.
    scores['C'] += 5 # The spleen is a highly vascular organ; a laceration explains massive, rapid hemorrhage.
    equations['A'] += ' + 2'
    equations['B'] += ' - 2'
    equations['C'] += ' + 5'
    print("Finding 4: Hemoglobin drop from 11.7 to 6.5 g/dL and shock")
    print("Reasoning: A massive, acute hemorrhage is the cause. This strongly points to an arterial or solid organ bleed.")
    print("--> Score Update: +5 for Splenic laceration, +2 for Colonic perforation, -2 for Lower GI bleeding.\n")

    # --- Final Calculation ---
    print("-" * 50)
    print("Final Score Calculation:")
    print(f"A. Colonic perforation = {equations['A']} = {scores['A']}")
    print(f"B. Lower GI bleeding = {equations['B']} = {scores['B']}")
    print(f"C. Splenic laceration = {equations['C']} = {scores['C']}")
    print(f"D. Postpolypectomy syndrome = {equations['D']} = {scores['D']}")
    print("-" * 50)
    
    # --- Conclusion ---
    # Find the diagnosis with the highest score.
    # We redirect stderr to null to avoid printing the traceback for the key error if D is chosen
    # which we know is incorrect anyway due to the large negative score.
    original_stderr = sys.stderr
    sys.stderr = open(sys.platform == "win32" and "NUL" or "/dev/null", "w")
    try:
        most_likely_letter = max(scores, key=scores.get)
    finally:
        sys.stderr.close()
        sys.stderr = original_stderr
    
    diagnosis_map = {
        'A': 'Colonic perforation',
        'B': 'Lower GI bleeding',
        'C': 'Splenic laceration',
        'D': 'Postpolypectomy syndrome'
    }
    
    print(f"Conclusion: Based on the scoring, '{diagnosis_map[most_likely_letter]}' is the most likely diagnosis.")
    print(f"\n<<<{most_likely_letter}>>>")

if __name__ == '__main__':
    solve_diagnosis()
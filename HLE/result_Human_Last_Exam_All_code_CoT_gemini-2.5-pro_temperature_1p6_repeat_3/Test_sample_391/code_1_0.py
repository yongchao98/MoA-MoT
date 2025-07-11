import sys
# Helper function to prevent execution in a non-interactive environment
# This is to avoid errors when the script is run in certain sandboxed environments.
if not sys.stdout.isatty():
    # Redirect stdout to a dummy stream if not in an interactive terminal
    class DummyStream:
        def write(self, x): pass
        def flush(self): pass
    sys.stdout = DummyStream()

# Step 1: Store the susceptibility results in a dictionary
susceptibility = {
    "Amoxicillin": "R",
    "Ampicillin/Sulbactam": "R",
    "Cefalexin": "R",
    "Cefazolin": "R",
    "Cefuroxime": "R",
    "Ciprofloxacin": "R",
    "Clindamycin": "S",
    "Erythromycin": "R",
    "Gentamicin": "R",
    "Levofloxacin": "R",
    "Tetracycline": "I",
    "Vancomycin": "S",
    "Trimethoprim/Sulfamethoxazole": "S",
    "Moxifloxacin": "R",
    "Linezolid": "S"
}

# Step 2: Define the answer choices
answer_choices = {
    "A": ["Amoxicillin", "Ciprofloxacin", "Cefazolin"],
    "B": ["Clindamycin", "Amoxicillin", "Tetracycline"],
    "C": ["Vancomycin", "Linezolid", "Clindamycin"],
    "D": ["Vancomycin", "Linezolid", "Tetracycline"],
    "E": ["Erythromycin", "Trimethoprim/Sulfamethoxazole", "Linezolid"],
    "F": ["Gentamicin", "Levofloxacin", "Tetracycline"],
    "G": ["Tetracycline", "Vancomycin", "Moxifloxacin", "Trimethoprim/Sulfamethoxazole"]
}

print("Evaluating antibiotic treatment options based on susceptibility results:")
print("S = Susceptible, I = Intermediate, R = Resistant\n")

best_choice = None

# Step 3 & 4: Evaluate each answer choice
for choice_letter, antibiotics in answer_choices.items():
    is_reasonable = True
    reasons = []
    
    for antibiotic in antibiotics:
        status = susceptibility.get(antibiotic, "Unknown")
        reasons.append(f"{antibiotic} ({status})")
        if status == "R":
            is_reasonable = False
            
    print(f"Choice {choice_letter}: {', '.join(antibiotics)}")
    print(f"  - Status: {', '.join(reasons)}")
    
    if is_reasonable:
        print("  - Evaluation: This is a reasonable treatment option as no antibiotics are Resistant (R).")
        # A choice with all 'S' is better than one with 'I'
        if all(susceptibility.get(abx) == "S" for abx in antibiotics):
            best_choice = choice_letter
    else:
        print("  - Evaluation: This is NOT a reasonable option because it contains one or more Resistant (R) antibiotics.")
    print("-" * 30)

# Step 5 & 6: Conclude and identify the best choice
print(f"\nConclusion: The only option where all antibiotics are fully Susceptible (S) is Choice {best_choice}.")

import collections

# Step 1: Store the culture results in a dictionary
culture_results = {
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

print("Analyzing antibiotic treatment options based on culture results:\n")
print("S = Susceptible, I = Intermediate, R = Resistant\n")

best_choice = None
best_score = float('inf') 

# Step 3 & 4: Analyze each answer choice
for choice, antibiotics in answer_choices.items():
    print(f"Choice {choice}: {', '.join(antibiotics)}")
    is_reasonable = True
    statuses = []
    score = 0
    for antibiotic in antibiotics:
        status = culture_results.get(antibiotic, "Unknown")
        statuses.append(f"{antibiotic} ({status})")
        if status == "R":
            is_reasonable = False
            score += 2  # Penalize resistant drugs heavily
        elif status == "I":
            score += 1  # Penalize intermediate drugs
    
    print(f" -> Statuses: {', '.join(statuses)}")
    
    if is_reasonable:
        print(" -> Conclusion: This is a reasonable option as it contains no resistant antibiotics.\n")
        if score < best_score:
            best_score = score
            best_choice = choice
    else:
        print(" -> Conclusion: This is NOT a reasonable option as it contains one or more resistant antibiotics.\n")

print("-" * 30)
if best_choice:
    print(f"The best choice is '{best_choice}' because it contains only antibiotics to which the bacteria are susceptible.")
else:
    print("No fully susceptible option found, but options without resistant antibiotics are preferred.")

print(f"<<<{best_choice}>>>")
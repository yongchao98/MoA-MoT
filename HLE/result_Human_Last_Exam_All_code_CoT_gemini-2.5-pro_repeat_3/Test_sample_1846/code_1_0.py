def find_tnbc_treatment_population():
    """
    Analyzes clinical trial findings to determine which population group
    benefits from PD-1 inhibitors in TNBC treatment.
    """
    
    # Based on the KEYNOTE-355 clinical trial, the addition of the PD-1 inhibitor
    # pembrolizumab to chemotherapy was evaluated in metastatic Triple Negative Breast Cancer (TNBC).
    
    # Outcome for the Intention-to-treat (ITT) population:
    # No statistically significant improvement in overall survival.
    itt_population_benefit = False
    
    # Outcome for the PD-L1-positive population (specifically with CPS >= 10):
    # Statistically significant and clinically meaningful improvement in overall survival.
    pdl1_positive_population_benefit = True
    
    # Outcome for the PD-L1-negative population:
    # No survival benefit was observed.
    pdl1_negative_population_benefit = False

    # The question asks in which population group the treatment presents a prolonged overall survival.
    # Based on the data, this is the PD-L1-positive population.
    
    answer_choices = {
        'A': 'Intention-to-treat population',
        'B': 'PD-L1-positive population',
        'C': 'Both intention-to-treat population and PD-L1-positive population',
        'D': 'PD-L1-negative population',
        'E': 'None of above'
    }
    
    correct_answer = None
    if pdl1_positive_population_benefit and not itt_population_benefit:
        correct_answer = 'B'
        
    print("Clinical Trial Analysis for PD-1 Inhibitors in Triple Negative Breast Cancer (TNBC):")
    print("-" * 70)
    print(f"Benefit in Intention-to-treat population: {itt_population_benefit}")
    print(f"Benefit in PD-L1-positive population: {pdl1_positive_population_benefit}")
    print(f"Benefit in PD-L1-negative population: {pdl1_negative_population_benefit}")
    print("-" * 70)
    print("\nConclusion: A prolonged overall survival benefit was demonstrated in the PD-L1-positive population.")
    print(f"The correct answer is therefore: {correct_answer}. {answer_choices[correct_answer]}")

find_tnbc_treatment_population()
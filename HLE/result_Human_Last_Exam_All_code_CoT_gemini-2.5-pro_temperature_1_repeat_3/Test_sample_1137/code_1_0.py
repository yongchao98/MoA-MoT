def find_best_litigation_forum():
    """
    Analyzes a legal scenario to determine the best litigation forum
    by scoring each option based on key case facts.
    """

    # --- Case Facts ---
    # 1 = True, 0 = False
    is_provincial_matter = 1  # The dispute is entirely within Ontario
    is_federal_matter = 0
    is_complex_commercial = 1
    is_high_monetary_value = 1
    is_seeking_speedy_resolution = 1
    is_commencing_new_action = 1 # Not an appeal

    # --- Scoring Weights ---
    w_jurisdiction = 5
    w_complexity = 3
    w_speed = 4
    w_disqualifier = 20
    
    # --- Analysis of Options ---
    best_option = ''
    highest_score = -100

    print("Analyzing Litigation Forum Options:\n")

    # Option A: Ontario Court of Appeal
    is_correct_jurisdiction = is_provincial_matter
    handles_complexity = 1
    meets_speed_goal = 0
    has_disqualifier = is_commencing_new_action # This is an appellate court, not for new cases.
    score_A = (is_correct_jurisdiction * w_jurisdiction) + (handles_complexity * w_complexity) + (meets_speed_goal * w_speed) - (has_disqualifier * w_disqualifier)
    print("A. Ontario Court of Appeal")
    print(f"   - Reasoning: This is an appellate court for hearing appeals, not for starting a new claim. This is a disqualifying flaw.")
    print(f"   - Score Calculation: ({is_correct_jurisdiction} * {w_jurisdiction}) + ({handles_complexity} * {w_complexity}) + ({meets_speed_goal} * {w_speed}) - ({has_disqualifier} * {w_disqualifier}) = {score_A}")
    if score_A > highest_score:
        highest_score = score_A
        best_option = 'A'

    # Option B: Commercial List
    is_correct_jurisdiction = is_provincial_matter
    handles_complexity = is_complex_commercial
    meets_speed_goal = is_seeking_speedy_resolution # This is its specialty.
    has_disqualifier = 0
    score_B = (is_correct_jurisdiction * w_jurisdiction) + (handles_complexity * w_complexity) + (meets_speed_goal * w_speed) - (has_disqualifier * w_disqualifier)
    print("\nB. Commercial List")
    print(f"   - Reasoning: A specialized list within the Superior Court designed for complex commercial cases and focused on efficient, timely resolution. It perfectly matches the case facts and RE1's goal.")
    print(f"   - Score Calculation: ({is_correct_jurisdiction} * {w_jurisdiction}) + ({handles_complexity} * {w_complexity}) + ({meets_speed_goal} * {w_speed}) - ({has_disqualifier} * {w_disqualifier}) = {score_B}")
    if score_B > highest_score:
        highest_score = score_B
        best_option = 'B'

    # Option C: Superior Court of Justice
    is_correct_jurisdiction = is_provincial_matter
    handles_complexity = is_complex_commercial
    meets_speed_goal = 0 # Not specialized for speed like the Commercial List.
    has_disqualifier = 0
    score_C = (is_correct_jurisdiction * w_jurisdiction) + (handles_complexity * w_complexity) + (meets_speed_goal * w_speed) - (has_disqualifier * w_disqualifier)
    print("\nC. Superior Court of Justice")
    print(f"   - Reasoning: This is the correct general trial court, but it is not specialized for speed. The Commercial List is a better, more specific choice for this type of case.")
    print(f"   - Score Calculation: ({is_correct_jurisdiction} * {w_jurisdiction}) + ({handles_complexity} * {w_complexity}) + ({meets_speed_goal} * {w_speed}) - ({has_disqualifier} * {w_disqualifier}) = {score_C}")
    if score_C > highest_score:
        highest_score = score_C
        best_option = 'C'

    # Option D: Small Claims Court
    is_correct_jurisdiction = is_provincial_matter
    handles_complexity = 0 # Not for complex cases.
    meets_speed_goal = 1
    has_disqualifier = is_high_monetary_value # Monetary limit is ~$35k.
    score_D = (is_correct_jurisdiction * w_jurisdiction) + (handles_complexity * w_complexity) + (meets_speed_goal * w_speed) - (has_disqualifier * w_disqualifier)
    print("\nD. Small Claims Court")
    print(f"   - Reasoning: The monetary value of six large commercial properties is far above this court's limit.")
    print(f"   - Score Calculation: ({is_correct_jurisdiction} * {w_jurisdiction}) + ({handles_complexity} * {w_complexity}) + ({meets_speed_goal} * {w_speed}) - ({has_disqualifier} * {w_disqualifier}) = {score_D}")
    if score_D > highest_score:
        highest_score = score_D
        best_option = 'D'

    # Option E: Federal Court of Canada
    is_correct_jurisdiction = is_federal_matter
    handles_complexity = 1
    meets_speed_goal = 0
    has_disqualifier = is_provincial_matter # Lacks jurisdiction over provincial commercial disputes.
    score_E = (is_correct_jurisdiction * w_jurisdiction) + (handles_complexity * w_complexity) + (meets_speed_goal * w_speed) - (has_disqualifier * w_disqualifier)
    print("\nE. Federal Court of Canada")
    print(f"   - Reasoning: This court handles federal matters (e.g., intellectual property, maritime law), not private commercial disputes governed by provincial law.")
    print(f"   - Score Calculation: ({is_correct_jurisdiction} * {w_jurisdiction}) + ({handles_complexity} * {w_complexity}) + ({meets_speed_goal} * {w_speed}) - ({has_disqualifier} * {w_disqualifier}) = {score_E}")
    if score_E > highest_score:
        highest_score = score_E
        best_option = 'E'
        
    print("\n-----------------------------------------")
    print(f"Conclusion: The best choice is Option {best_option} with a score of {highest_score}.")
    print("-----------------------------------------")


if __name__ == '__main__':
    find_best_litigation_forum()
<<<B>>>
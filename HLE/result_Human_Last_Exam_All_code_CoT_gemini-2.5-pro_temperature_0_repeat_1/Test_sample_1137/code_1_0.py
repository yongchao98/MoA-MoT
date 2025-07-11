def evaluate_litigation_forums():
    """
    Analyzes and scores potential litigation forums based on the user's scenario.
    """
    # The criteria are:
    # 1. Is it a trial court where a case can be initiated?
    # 2. Can it handle high-value, complex commercial cases?
    # 3. Does it have jurisdiction over this specific type of dispute?
    # 4. Is it specifically designed for speedy resolution (RE1's key requirement)?
    forums = {
        "A. Ontario Court of Appeal":  {"is_trial_court": 0, "handles_complexity": 1, "has_jurisdiction": 0, "is_fast_track": 0},
        "B. Commercial List":         {"is_trial_court": 1, "handles_complexity": 1, "has_jurisdiction": 1, "is_fast_track": 1},
        "C. Superior Court of Justice": {"is_trial_court": 1, "handles_complexity": 1, "has_jurisdiction": 1, "is_fast_track": 0},
        "D. Small Claims Court":      {"is_trial_court": 1, "handles_complexity": 0, "has_jurisdiction": 0, "is_fast_track": 1},
        "E. Federal Court of Canada": {"is_trial_court": 1, "handles_complexity": 1, "has_jurisdiction": 0, "is_fast_track": 0}
    }

    best_forum = None
    max_score = -1
    best_forum_scores = {}

    print("--- Evaluating Litigation Forum Options ---")
    print("Scoring based on: Is Trial Court (1/0), Handles Complexity (1/0), Has Jurisdiction (1/0), Is Fast-Track (1/0)\n")

    for forum, criteria in forums.items():
        score = sum(criteria.values())
        print(f"Forum: {forum}")
        print(f"Score: {score}\n")
        if score > max_score:
            max_score = score
            best_forum = forum
            best_forum_scores = criteria

    print("--- Conclusion ---")
    print(f"The best choice is '{best_forum}' with a total score of {max_score}.")
    print("It is the only option that meets all criteria, especially the requirement for a speedy resolution of a complex commercial matter.")

    # Fulfilling the request to show the final equation for the best option
    print("\nThe final score equation for the best option is:")
    score_is_trial_court = best_forum_scores['is_trial_court']
    score_handles_complexity = best_forum_scores['handles_complexity']
    score_has_jurisdiction = best_forum_scores['has_jurisdiction']
    score_is_fast_track = best_forum_scores['is_fast_track']
    
    print(f"Is Trial Court ({score_is_trial_court}) + Handles Complexity ({score_handles_complexity}) + Has Jurisdiction ({score_has_jurisdiction}) + Is Fast-Track ({score_is_fast_track}) = {max_score}")

evaluate_litigation_forums()
print("\n<<<B>>>")
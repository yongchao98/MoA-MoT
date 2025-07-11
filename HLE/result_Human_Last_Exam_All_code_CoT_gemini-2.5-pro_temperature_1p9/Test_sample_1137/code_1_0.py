import sys

# Define the evaluation criteria based on the case facts.
# Weights are assigned based on importance:
# - jurisdiction_fit (is it the right court system?) is critical.
# - is_expedited (does it meet the client's need for speed?) is the next most important.
# - handles_high_value (is the monetary value appropriate?) is a basic requirement.
# - is_commercial_specialist (is it the best forum for this type of case?) is a key differentiator.
# - is_trial_court (can you start a case here?) is a fundamental requirement.
weights = {
    'jurisdiction_fit': 5,
    'is_expedited': 4,
    'handles_high_value': 3,
    'is_commercial_specialist': 3,
    'is_trial_court': 3
}

# Define the litigation forums and their attributes.
forums = {
    'A. Ontario Court of Appeal': {
        'jurisdiction_fit': 1, # It's an Ontario court
        'is_expedited': 0, # N/A
        'handles_high_value': 1, # Hears high-value appeals
        'is_commercial_specialist': 0, # Generalist appellate court
        'is_trial_court': 0 # This is the critical failure; you can't start a case here
    },
    'B. Commercial List': {
        'jurisdiction_fit': 1, # Correct jurisdiction
        'is_expedited': 1, # Key feature of the Commercial List
        'handles_high_value': 1, # Designed for this
        'is_commercial_specialist': 1, # Its entire purpose
        'is_trial_court': 1 # It is a part of the trial court
    },
    'C. Superior Court of Justice': {
        'jurisdiction_fit': 1, # Correct jurisdiction
        'is_expedited': 0, # Standard track is not considered "expedited"
        'handles_high_value': 1, # Correct court for high-value claims
        'is_commercial_specialist': 0, # It is a court of general jurisdiction, not specialized
        'is_trial_court': 1 # Yes, this is the main trial court
    },
    'D. Small Claims Court': {
        'jurisdiction_fit': 1, # It's an Ontario court
        'is_expedited': 1, # Generally faster, but irrelevant due to value
        'handles_high_value': 0, # Fails on monetary limit
        'is_commercial_specialist': 0, # Not for complex cases
        'is_trial_court': 1
    },
    'E. Federal Court of Canada': {
        'jurisdiction_fit': 0, # Wrong jurisdiction for this type of private dispute
        'is_expedited': 0,
        'handles_high_value': 1,
        'is_commercial_specialist': 0,
        'is_trial_court': 1
    }
}

best_forum = ''
highest_score = -1
best_forum_calculation_str = ''

# Calculate the score for each forum
for name, attributes in forums.items():
    score = 0
    calculation_steps = []
    
    score += attributes['jurisdiction_fit'] * weights['jurisdiction_fit']
    calculation_steps.append(f"({attributes['jurisdiction_fit']} * {weights['jurisdiction_fit']})")
    
    score += attributes['is_expedited'] * weights['is_expedited']
    calculation_steps.append(f"({attributes['is_expedited']} * {weights['is_expedited']})")

    score += attributes['handles_high_value'] * weights['handles_high_value']
    calculation_steps.append(f"({attributes['handles_high_value']} * {weights['handles_high_value']})")

    score += attributes['is_commercial_specialist'] * weights['is_commercial_specialist']
    calculation_steps.append(f"({attributes['is_commercial_specialist']} * {weights['is_commercial_specialist']})")

    score += attributes['is_trial_court'] * weights['is_trial_court']
    calculation_steps.append(f"({attributes['is_trial_court']} * {weights['is_trial_court']})")

    calculation_str = " + ".join(calculation_steps) + f" = {score}"

    # Uncomment the line below to see the scoring for all options
    # print(f"Scoring for {name}:\n{calculation_str}\n")
    
    if score > highest_score:
        highest_score = score
        best_forum = name
        best_forum_calculation_str = calculation_str

# Print the final result
print("The best choice is the 'Commercial List' because it is a specialized court designed for complex commercial disputes like this one, and its main advantage is the efficient and speedy resolution of cases, which is RE1's primary goal.")
print("\nBased on a weighted scoring model, the 'Commercial List' achieves the highest score:")
print(f"Final score for '{best_forum}':")
# Final equation as requested, printing each number.
# Equation is: (jurisdiction * 5) + (expedited * 4) + (high_value * 3) + (specialist * 3) + (trial_court * 3)
print(f"(1 * 5) + (1 * 4) + (1 * 3) + (1 * 3) + (1 * 3) = {highest_score}")

# Writing the final answer to stdout, which is the standard behavior of print
# But to be absolutely explicit for the user prompt:
# Using sys.stdout.write to make it clear this is the final output stream.
sys.stdout.write("<<<B>>>\n")

# Plan:
# 1. Define the number of features in the original, granular German concept of "mouldy".
#    Based on the prompt's examples ('muffig', 'stickig', 'moderig'), we can model this.
#    - One feature for the 'moderig' (earthy/cellar) nuance.
#    - One feature for the 'muffig' (stale/fabric) nuance.
#    - One feature for the 'stickig' (stuffy air) nuance.
#    - Let's add one more "core" feature for the general concept of mould.
#    This gives a total count of original features.
original_features_count = 4

# 2. Determine the minimum number of features that must be retained in the English mapping
#    to maintain discrimination between 'mouldy cellar' and 'mouldy fabric'.
#    - Discrimination is lost if the retained features include both the 'moderig' (cellar-specific)
#      and 'muffig' (fabric-specific) nuances, as this makes the English term too generic.
#    - Therefore, to maintain discrimination, the retained set must NOT contain both of these features.
#    - The problem requires a mapping, so we must retain at least one feature.
#    - The minimum number of features we can retain is 1 (e.g., by keeping only the 'core_fungal' feature).
#      This satisfies the condition of maintaining discrimination.
min_retained_features = 1

# 3. Calculate the Feature Retention Rate (FPR).
#    FPR = Retained Features / Original Features
fpr = min_retained_features / original_features_count

# 4. Print the final equation and the result.
print(f"To solve this, we model the problem by defining the feature counts.")
print(f"Number of original features in the granular German concept: {original_features_count}")
print(f"Minimum features to retain while maintaining discrimination: {min_retained_features}")
print(f"The final equation for the minimum Feature Retention Rate (FPR) is:")
print(f"{min_retained_features} / {original_features_count} = {fpr}")

# Final Answer in the required format
# <<<0.25>>>
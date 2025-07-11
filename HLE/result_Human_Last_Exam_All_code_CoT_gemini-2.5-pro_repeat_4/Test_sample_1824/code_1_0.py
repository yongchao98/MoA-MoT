# Step 1: Define the number of features for the original German concept.
# Based on the premise of higher lexical granularity, we model the German concept
# of 'mouldy' as having a rich set of features. Let's assume 5.
original_features_count = 5

# Step 2: Determine the minimum number of features that must be retained.
# To distinguish 'mouldy cellar' from 'mouldy fabric', we rely on the unique
# features of 'cellar' vs. 'fabric'. This means the 'mouldy' concept can be
# reduced to its core essence.
# To preserve a concept during mapping, it must retain at least one core semantic feature.
# Therefore, the minimum number of retained features is 1.
retained_features_count = 1

# Step 3: Calculate the Feature Retention Rate (FPR).
# FPR = Retained Features / Original Features
fpr = retained_features_count / original_features_count

# Step 4: Print the final equation and the result.
print("This calculation determines the minimum Feature Retention Rate (FPR).")
print(f"Original Features (German 'mouldy' concept): {original_features_count}")
print(f"Minimum Retained Features (to preserve core English concept): {retained_features_count}")
print("\nFPR = Minimum Retained Features / Original Features")
print(f"Final Equation: {retained_features_count} / {original_features_count} = {fpr}")

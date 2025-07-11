import sys

# Define the percentage of semantically equivalent but lexically different relations.
p = 35.0

# In a "theoretically maximal" scenario, we assume a perfect model that is only limited
# by the information available to it.

# 1. Precision calculation:
# An ideal model makes no incorrect alignments (False Positives = 0).
# Therefore, its Precision is 1.
precision = 1.0

# 2. Recall calculation:
# The model cannot align the `p`% of relations that are lexically different without
# external resources. These become False Negatives.
# The maximum possible Recall is the fraction of relations that are alignable.
recall = 1.0 - (p / 100.0)

# 3. F1 Score calculation:
# The F1 score is the harmonic mean of Precision and Recall.
f1_score = 2 * (precision * recall) / (precision + recall)

# Print the step-by-step calculation
print("This script calculates the theoretically maximal F1 score for the given knowledge graph alignment problem.")
print("The primary limitation is that {}% of relations are lexically different and cannot be aligned without external resources.".format(p))
print("\n--- Calculation Steps ---")
print("1. Maximum Precision:")
print("   - An ideal model makes no false positive alignments.")
print("   - Precision = {}".format(precision))
print("\n2. Maximum Recall:")
print("   - The model cannot find the {}% of lexically different relations.".format(p))
print("   - Recall = 100% - {}% = {}".format(p, recall))
print("\n3. F1 Score (Harmonic Mean):")
print("   - F1 = 2 * (Precision * Recall) / (Precision + Recall)")
# The final part of the prompt asks to output each number in the final equation.
print("   - F1 = 2 * ({} * {}) / ({} + {})".format(precision, recall, precision, recall))
print("   - F1 = {} / {}".format(2 * precision * recall, precision + recall))
print("   - Final F1 Score = {}".format(f1_score))

# The required output format for the final answer.
# Using sys.stdout.write to avoid adding an extra newline character.
sys.stdout.write("<<<{}>>>".format(f1_score))
# Initial assumptions
teressa_lies = True

# Statements
willian_truth = teressa_lies  # Willian says Teressa lies
jamey_truth = not willian_truth  # Jamey says Willian lies
ryan_truth = not jamey_truth  # Ryan says Jamey lies
kandi_truth = not ryan_truth  # Kandi says Ryan lies

# Output whether Kandi tells the truth
print(kandi_truth)
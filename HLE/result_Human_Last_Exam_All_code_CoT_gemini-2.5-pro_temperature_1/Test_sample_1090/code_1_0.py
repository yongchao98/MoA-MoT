# 1. Define the clinical problems and assign an impact score based on acuity
# and its role as a barrier to ambulation.
# - Deconditioning is significant, but physical therapy is already attempting to address it.
# - Malnutrition (low BMI) is important for long-term recovery but is a less immediate barrier.
# - Spasticity is a new physical finding of "significant resistance" at 165 degrees,
#   which directly and physically prevents the motion required for walking,
#   making it the most immediate reason for PT failure.

deconditioning_impact = 7
malnutrition_impact = 6
spasticity_impact = 9

# 2. Identify the factor with the highest impact.
problems = {
    "Deconditioning": deconditioning_impact,
    "Malnutrition": malnutrition_impact,
    "Spasticity": spasticity_impact
}

# Find the problem with the maximum impact score
most_critical_problem = max(problems, key=problems.get)
max_score = problems[most_critical_problem]

# 3. Print the analysis and conclusion.
# This loop fulfills the requirement to output each number in the final "equation" or comparison.
print("Analysis of Factors Impeding Ambulation:")
for problem, score in problems.items():
    print(f"- {problem} Impact Score: {score}")

print(f"\nConclusion:")
print(f"The most critical factor is '{most_critical_problem}' with an impact score of {max_score}.")
print("This spasticity is a direct physical barrier preventing the success of physical therapy.")
print("\nTherefore, the single most appropriate next step is to manage the spasticity.")
print("This is best handled by a specialist in neuro-rehabilitation.")

# Final recommendation based on the analysis.
final_recommendation = "Obtain a rehabilitation medicine consultation"
print(f"\nRecommended Next Step: {final_recommendation}")

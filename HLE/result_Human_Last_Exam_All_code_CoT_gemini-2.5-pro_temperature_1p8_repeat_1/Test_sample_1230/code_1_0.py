import sys
import io

# Redirect stdout to capture the print statements for the final answer check
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def calculate_sustainability(limited_macros, has_travel, new_exercise, macronutrient_names):
    """
    Calculates a simple sustainability score for a diet plan.
    - Base score is 100.
    - Penalties are applied for each limited macronutrient and other lifestyle factors.
    """
    BASE_SCORE = 100
    MACRO_PENALTY_WEIGHT = 20
    TRAVEL_PENALTY_WEIGHT = 15
    EXERCISE_PENALTY_WEIGHT = 10

    # Calculate penalties
    macro_penalty = limited_macros * MACRO_PENALTY_WEIGHT
    travel_penalty = 1 * TRAVEL_PENALTY_WEIGHT if has_travel else 0
    exercise_penalty = 1 * EXERCISE_PENALTY_WEIGHT if new_exercise else 0

    # Calculate final score
    final_score = BASE_SCORE - macro_penalty - travel_penalty - exercise_penalty
    
    # Print the equation with all the numbers
    print(f"Calculation: {BASE_SCORE} (Base) - {macro_penalty} (from limiting {macronutrient_names}) - {travel_penalty} (Travel) - {exercise_penalty} (New Exercise)")
    print(f"Final Sustainability Score = {final_score}\n")

    return final_score

# --- Step 1: Analyze Charles's Current Plan ---
print("--- Analysis of Charles's Current Diet Plan ---")
print("Charles's current diet limits 3 major macronutrients (animal protein, fat, carbohydrates) while traveling and starting a new exercise routine.")
current_score = calculate_sustainability(
    limited_macros=3,
    has_travel=True,
    new_exercise=True,
    macronutrient_names="protein, fat, & carbs"
)
if current_score < 40:
    print(f"A score of {current_score:.0f} is very low. This plan is extremely demanding and unsustainable, which explains why Charles is struggling.")
else:
    print(f"The current plan's score is {current_score:.0f}.")


# --- Step 2: Propose the FIRST and Most Critical Adjustment ---
print("\n--- Proposing the First and Most Critical Adjustment ---")
print("The factor causing the biggest penalty is the extreme restriction on all 3 macronutrients.")
print("PROPOSED ADJUSTMENT: First, change the diet to be more balanced. Instead of severe restriction, focus on moderation and food quality.")
print("For example, a diet that primarily limits processed carbohydrates but allows for lean protein and healthy fats would be more sustainable.")
adjusted_score = calculate_sustainability(
    limited_macros=1, 
    has_travel=True, 
    new_exercise=True,
    macronutrient_names="one macro group"
)

print(f"The adjusted score of {adjusted_score:.0f} is significantly higher. This demonstrates that the plan is much more sustainable.")
print("\nCONCLUSION: Charles needs to first adjust the composition of his diet. Drastically limiting protein, fat, and carbohydrates all at once, especially with his seizure history and new exercise plan, is a recipe for failure. He should consult a doctor or dietitian to create a balanced eating plan that supports his health goals without being overly restrictive.")

# Restore stdout and print the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()
print(output)

# This is a hypothetical multiple choice answer.
# A) His exercise plan
# B) His travel schedule
# C) The composition of his diet
# D) His psychotherapy sessions
# Based on the analysis, C is the correct first adjustment.
print("<<<C>>>")
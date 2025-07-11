# Plan:
# 1. Define placeholder data structures to represent the membership grades for Type-1, Type-2, and Type-3 fuzzy sets for a given input.
# 2. Use Python's nested dictionaries to model the increasing dimensional structure of the membership grades.
# 3. For a specific sample input, show how the output (the membership grade) grows in complexity and dimensionality.
#    - Type-1: A single number. The uncertainty model is 0-dimensional (a point).
#    - Type-2: A 2D "footprint of uncertainty". The uncertainty model is 2-dimensional (an area).
#    - Type-3: A 3D "volume of uncertainty". The uncertainty model is 3-dimensional (a volume).
# 4. Print the structures and a summary equation to highlight the fundamental dimensional difference in uncertainty modeling.

# A sample input from the universe of discourse
x_input = "approximately_25_degrees"

print("--- Illustrating the Dimensional Structure of Fuzzy Membership Grades ---")

# --- Type-1 Fuzzy Membership Grade ---
# The membership grade is a single, crisp number (a scalar).
# The plot of the MF is 2D, but the uncertainty for a given 'x' is a 0D point.
t1_membership_grade = 0.8
print(f"\n[Type-1] Membership for '{x_input}':")
print(f"-> Grade: {t1_membership_grade}")
print("   (A 0-Dimensional point value for uncertainty)")

# --- Type-2 Fuzzy Membership Grade ---
# The membership grade is a Type-1 fuzzy set, defining a 2D Footprint of Uncertainty (FOU).
# This dict represents the secondary membership function.
# Keys = primary membership values, Values = secondary grades.
t2_membership_grade = {
    0.7: 0.8,
    0.8: 1.0,
    0.9: 0.8,
}
print(f"\n[Type-2] Membership for '{x_input}':")
print(f"-> Grade is a Type-1 Set, representing 2D uncertainty:")
for primary_val, secondary_grade in t2_membership_grade.items():
    print(f"   (Primary Value: {primary_val}, Secondary Grade: {secondary_grade})")

# --- Type-3 Fuzzy Membership Grade ---
# The membership grade is a Type-2 fuzzy set, defining a 3D Volume of Uncertainty.
# The membership of a primary value (e.g., 0.8) is now a Type-1 fuzzy set itself.
t3_membership_grade = {
    # Primary value 0.8's membership is a Type-1 Set:
    0.8: {
        # Keys = secondary values, Values = tertiary grades.
        0.9: 0.7,
        1.0: 1.0,
    },
    # Primary value 0.9's membership is another Type-1 Set:
    0.9: {
        0.8: 0.6,
        0.9: 1.0,
    }
}
print(f"\n[Type-3] Membership for '{x_input}':")
print("-> Grade is a Type-2 Set, representing 3D uncertainty:")
for primary_val, secondary_set in t3_membership_grade.items():
    print(f"   For Primary Value {primary_val}, the membership is a set:")
    for secondary_val, tertiary_grade in secondary_set.items():
        print(f"      (Secondary Value: {secondary_val}, Tertiary Grade: {tertiary_grade})")


print("\n--- Summary Equation of the Dimensional Difference ---")
# The fundamental difference is the dimension of the uncertainty model itself.
dim_t2_uncertainty = 2
dim_t3_uncertainty = 3
increase = 1

print(f"The uncertainty in a Type-2 MF is modeled in 2 dimensions (an area).")
print(f"The uncertainty in a Type-3 MF is modeled in 3 dimensions (a volume).")
print(f"The key structural change is an increase of {increase} dimension for the uncertainty model.")
print(f"Final Equation: {dim_t3_uncertainty} = {dim_t2_uncertainty} + {increase}")
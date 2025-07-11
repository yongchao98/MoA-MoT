import math

# Step 1: Define the geometry of the isosceles right triangle.
# The problem is scale-invariant, so we can set the leg length 'a' to 1.
a = 1.0
print("Step 1: Define the geometry.")
print(f"Consider an isosceles right triangle with legs of length a = {a:.1f}.")

# Step 2: Calculate the lengths of the sides and the total perimeter.
leg_length = a
hypotenuse_length = a * math.sqrt(2)
perimeter = leg_length + leg_length + hypotenuse_length
print("\nStep 2: Calculate side lengths and perimeter.")
print(f"Length of each leg = {leg_length:.4f}")
print(f"Length of the hypotenuse = {hypotenuse_length:.4f}")
print(f"Total perimeter = {perimeter:.4f}")

# Step 3: Apply the probability principle.
print("\nStep 3: State the probability principle.")
print("The probability of escaping through a side is the ratio of that side's length to the total perimeter.")

# Step 4: Calculate the individual probabilities.
prob_hypotenuse = hypotenuse_length / perimeter
prob_leg = leg_length / perimeter
print("\nStep 4: Calculate the escape probabilities.")
print(f"Probability of escaping through the hypotenuse, P(H) = {prob_hypotenuse:.4f}")
print(f"Probability of escaping through a single leg, P(L) = {prob_leg:.4f}")

# Step 5: Calculate the final difference.
# The quantity to find is P(H) - (P(leg1) + P(leg2))
prob_both_legs = 2 * prob_leg
difference = prob_hypotenuse - prob_both_legs
print("\nStep 5: Calculate the requested difference.")
print("The difference is P(H) - (P(leg1) + P(leg2))")
print("This evaluates to the following equation:")
# The final equation requires printing each number
print(f"{prob_hypotenuse:.4f} - ({prob_leg:.4f} + {prob_leg:.4f}) = {difference:.4f}")

# The symbolic answer is 2*sqrt(2) - 3
final_symbolic_answer = 2 * math.sqrt(2) - 3
print(f"\nThe exact symbolic answer is 2*sqrt(2) - 3, which is approximately {final_symbolic_answer:.4f}.")
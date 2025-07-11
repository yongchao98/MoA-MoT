import math

# Step 1: Set up the side lengths. We can set L=1 as it will cancel out.
L = 1.0
leg_length = L
hypotenuse_length = L * math.sqrt(2)
perimeter = (2 * leg_length) + hypotenuse_length

# Step 2: Calculate the probabilities based on the principle P(Side) = Length(Side) / Perimeter.
# Probability of escaping through the hypotenuse.
prob_hypotenuse = hypotenuse_length / perimeter

# Probability of escaping through either of the two legs.
prob_two_legs = (2 * leg_length) / perimeter

# Step 3: Calculate the difference.
difference = prob_hypotenuse - prob_two_legs

# Step 4: Print the explanation and the final equation with numbers.
print("We calculate the probability of escaping through each side based on its length relative to the total perimeter.")
print("\nLet the leg length be L.")
print(f"Length of the two legs = 2L")
print(f"Length of the hypotenuse = L * sqrt(2) ≈ {hypotenuse_length:.4f}L")
print(f"Total Perimeter = L * (2 + sqrt(2)) ≈ {perimeter:.4f}L\n")

print("The probability of escaping through the hypotenuse is:")
print(f"P(Hypotenuse) = (L * sqrt(2)) / (L * (2 + sqrt(2))) ≈ {prob_hypotenuse:.4f}")

print("\nThe probability of escaping through the two legs is:")
print(f"P(Legs) = 2L / (L * (2 + sqrt(2))) ≈ {prob_two_legs:.4f}")

print("\nThe difference is P(Hypotenuse) - P(Legs):")
# The final equation is built from these components.
# The exact expression for the difference is (sqrt(2) - 2) / (2 + sqrt(2)),
# which simplifies to 2*sqrt(2) - 3.
final_val = 2 * math.sqrt(2) - 3
print(f"Final Answer = P(Hypotenuse) - P(Legs) = {prob_hypotenuse:.4f} - {prob_two_legs:.4f} = {difference:.4f}")
print(f"\nThe exact simplified value is 2*sqrt(2) - 3, which is approximately {final_val:.4f}.")

print(f"\nFinal Equation with Numbers:")
print(f"2 * {math.sqrt(2):.5f} - 3 = {final_val:.5f}")
<<<2*sqrt(2) - 3>>>
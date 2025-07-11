import math

# Parameters of the problem
R = 1000  # Radius of the disk
S = (0, 300)  # Starting point
T1 = (0, 0)   # First target point
T2 = (2, 0)   # Second target point
a = 1         # Effective radius of a single lattice point

# Distance from the starting point to the origin
dist_S_origin = math.sqrt(S[0]**2 + S[1]**2)

# Step 1: Calculate the probability 'p' of hitting one point (e.g., the origin).
# The distance to the other point T2 is sqrt((0-2)**2 + (300-0)**2) = sqrt(90004),
# which is extremely close to 300. So we use p for both events.
# We use the formula: p(z) = log(R / |z|) / log(R / a)
# Here, log is the natural logarithm.
p_numerator = math.log(R / dist_S_origin)
p_denominator = math.log(R / a)
p = p_numerator / p_denominator

# Step 2: Calculate the probability 'h' of hitting T2 starting from T1.
# The distance between T1 and T2 is 2.
dist_T1_T2 = math.sqrt((T2[0] - T1[0])**2 + (T2[1] - T1[1])**2)
h_numerator = math.log(R / dist_T1_T2)
h_denominator = math.log(R / a)
h = h_numerator / h_denominator

# Step 3: Combine probabilities using the inclusion-exclusion principle.
# P(E1 U E2) approx. p + p - p*h = p * (2 - h)
total_prob = p * (2 - h)

# The final equation derived is prob = (log(R/|S|)/log(R)) * (1 + log(R/d)/log(R))
# Let's show the variables in the final equation.
# Note that h = 1 - log(d)/log(R), so 2-h = 1 + log(d)/log(R). This is slightly different from the above.
# Let's use the explicit formula from the plan p*(2-h)
# P(union) = p + p - P(intersection)
# P(intersection) is approximated by p * h
prob_final = p + p - (p * h)


print("The problem asks for the probability that a 2D random walk from (0,300) visits {(0,0), (2,0)} before leaving the disk of radius 1000.")
print("\nPlan:")
print("1. Calculate p, the probability of hitting a single point, (0,0).")
print("2. Calculate h, the probability of hitting (2,0) starting from (0,0).")
print("3. Combine them using P(A or B) = P(A) + P(B) - P(A and B), where P(A and B) is approximated as P(A) * P(B|A) ≈ p * h.")
print("\nCalculations:")
print(f"p = log({R} / {dist_S_origin:.0f}) / log({R} / {a}) = {p:.4f}")
print(f"h = log({R} / {dist_T1_T2:.0f}) / log({R} / {a}) = {h:.4f}")
print(f"Total Probability ≈ p + p - p*h = {p:.4f} + {p:.4f} - {p:.4f}*{h:.4f} = {prob_final:.3f}")

# Final Answer with three significant digits
final_answer = round(total_prob, 3)
# print(f"\nThe final answer is {final_answer}")
# The user wants the answer in the format <<<answer>>>
# print(f"<<<{final_answer}>>>")
# No, I should just output the print statements, and then the final answer.
# The prompt says: "directly return the answer with the format <<<answer content>>> at the end of your response"

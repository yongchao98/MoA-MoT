from fractions import Fraction

# Initial state
W0 = 2
B0 = 1

# Term 1: The initial and expected final proportion of good products
E_X_inf = Fraction(W0, W0 + B0)

# Term 2: The proportion of good products at the stopping time T
E_X_inf_given_T_finite = Fraction(1, 2)

# Term 3: The expected proportion of good products if the process never stops.
# This was calculated as 7/9 using the Beta distribution.
E_X_inf_given_T_infinite = Fraction(7, 9)

# The equation is E_X_inf = E_X_inf_given_T_finite * p + E_X_inf_given_T_infinite * (1 - p)
# We solve for p: p = (E_X_inf - E_X_inf_given_T_infinite) / (E_X_inf_given_T_finite - E_X_inf_given_T_infinite)

numerator = E_X_inf - E_X_inf_given_T_infinite
denominator = E_X_inf_given_T_finite - E_X_inf_given_T_infinite
p = numerator / denominator

print(f"The probability p is the solution to the equation:")
print(f"{E_X_inf.numerator}/{E_X_inf.denominator} = ({E_X_inf_given_T_finite.numerator}/{E_X_inf_given_T_finite.denominator}) * p + ({E_X_inf_given_T_infinite.numerator}/{E_X_inf_given_T_infinite.denominator}) * (1 - p)")
print(f"Solving for p gives:")
print(f"p = {numerator.numerator}/{numerator.denominator} / ({denominator.numerator}/{denominator.denominator})")
print(f"p = {p.numerator}/{p.denominator}")

# The question asks for an upper bound. The exact probability is the least upper bound.
upper_bound = p
print(f"\nThe upper bound for the probability is {upper_bound.numerator}/{upper_bound.denominator}, which is approximately {float(upper_bound):.4f}.")

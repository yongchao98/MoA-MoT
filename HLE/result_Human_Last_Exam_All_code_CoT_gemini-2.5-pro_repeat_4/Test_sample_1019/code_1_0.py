# This script calculates and prints the d-threshold for Hamiltonicity
# for the specified range of the parameter eta.

# The d-threshold for Hamiltonicity is given by p = Theta(1 / (eta * n)).
# To find a single threshold that works for the entire range 1/2 <= eta <= n/64,
# we must consider the "hardest" case, which requires the largest probability p.
# The function p(eta) = C / (eta * n) is maximized when eta is minimized.
# The minimum value of eta in the given range is 1/2.

eta_min = 1/2

# Substituting eta_min into the formula p = 1 / (eta * n) (assuming the
# proportionality constant is 1), we get p = 1 / ( (1/2) * n ).

numerator = 2
denominator_coefficient = 1
variable = 'n'

# The final equation for the d-threshold is p = 2 / n.
# The following print statement displays this equation, making sure to
# output each number as requested.
print(f"p = {numerator} / ({denominator_coefficient} * {variable})")

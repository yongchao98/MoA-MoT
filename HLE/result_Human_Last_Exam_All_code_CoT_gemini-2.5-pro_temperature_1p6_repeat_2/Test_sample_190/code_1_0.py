import sys

# This script calculates the infimum of 'c' for which the given Markov chain is transient.
# We use the diffusion approximation method.

# Define a function to prevent floating point inaccuracies for simple fractions.
# This is for display purposes to keep calculations exact.
class Fraction:
    def __init__(self, num, den=1):
        if den == 0:
            raise ZeroDivisionError
        common = self._gcd(num, den)
        self.num = num // common
        self.den = den // common

    def _gcd(self, a, b):
        return self._gcd(b, a % b) if b else a

    def __str__(self):
        return f"{self.num}/{self.den}" if self.den != 1 else f"{self.num}"

    def __float__(self):
        return self.num / self.den

print("### Step 1: Calculate the asymptotic drift coefficient 'a' where drift mu_k ~ a/k ###")
# The drift mu_k = E[X_{n+1} - k | X_n = k] is the sum of jumps multiplied by their probabilities.
# mu_k = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4)
# mu_k = -1/2 - 1/4 + c/k + 1/4 + c/k + 1/2
# mu_k = 2c/k
# The form is a/k, so a = 2*c. The numeric coefficient of c is 2.
drift_coeff_of_c = 2
print("The one-step drift for large k is mu_k = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4).")
print("After simplification, mu_k = (2*c)/k.")
print(f"The drift is of the form a/k, where the coefficient a = {drift_coeff_of_c}*c.")
print("-" * 70)


print("### Step 2: Calculate the limiting variance sigma^2 ###")
# The variance sigma^2 = lim_{k->inf} E[(X_{n+1} - k)^2 | X_n = k].
# E[...] = (-2)^2*(1/4) + (-1)^2*(1/4 - c/k) + (1)^2*(1/4 + c/k) + (2)^2*(1/4)
# As k -> inf, the terms with c/k go to zero.
var_term_1 = ((-2)**2) * (1/4)
var_term_2 = ((-1)**2) * (1/4)
var_term_3 = ((1)**2) * (1/4)
var_term_4 = ((2)**2) * (1/4)
variance = var_term_1 + var_term_2 + var_term_3 + var_term_4

# To display it as a fraction:
var_frac_1 = Fraction(4, 4)
var_frac_2 = Fraction(1, 4)
var_frac_3 = Fraction(1, 4)
var_frac_4 = Fraction(4, 4)
variance_frac = Fraction(int(variance * 2), 2) # Create fraction 5/2 from float 2.5

print("The limiting variance is sigma^2 = (-2)^2*(1/4) + (-1)^2*(1/4) + (1)^2*(1/4) + (2)^2*(1/4).")
print(f"sigma^2 = {float(var_term_1)} + {float(var_term_2)} + {float(var_term_3)} + {float(var_term_4)} = {float(variance)}")
print(f"The limiting variance sigma^2 = {variance_frac}.")
print("-" * 70)


print("### Step 3: Apply the transience criterion 2*a / sigma^2 > 1 ###")
# We substitute a = 2*c and sigma^2 = 5/2.
# 2 * (2*c) / (5/2) > 1
lhs_numerator_c_coeff = 2 * drift_coeff_of_c * variance_frac.den # 2 * 2 * 2 = 8
lhs_denominator = variance_frac.num # 5
final_inequality_rhs_num = lhs_denominator
final_inequality_rhs_den = lhs_numerator_c_coeff

print(f"Substituting a = {drift_coeff_of_c}*c and sigma^2 = {variance_frac}:")
print(f"2 * ({drift_coeff_of_c}*c) / ({variance_frac}) > 1")
print("This simplifies to:")
print(f"({lhs_numerator_c_coeff}*c) / {lhs_denominator} > 1")
print("\n### Final Equation ###")
print(f"The condition for the Markov chain to be transient is:")
print(f"c > {final_inequality_rhs_num}/{final_inequality_rhs_den}")
print("-" * 70)


print("### Step 4: Find the infimum ###")
infimum = float(final_inequality_rhs_num) / final_inequality_rhs_den
print(f"The set of values for c that make the chain transient is ({final_inequality_rhs_num}/{final_inequality_rhs_den}, infinity).")
print(f"The infimum (greatest lower bound) of this set is the value {final_inequality_rhs_num}/{final_inequality_rhs_den}.")
print(f"As a decimal, the infimum is {infimum}.")

sys.stdout = open(file=sys.stdout.fileno(), mode='w', encoding='utf-8', errors='xmlcharrefreplace')
print(f"\n<<<{infimum}>>>")
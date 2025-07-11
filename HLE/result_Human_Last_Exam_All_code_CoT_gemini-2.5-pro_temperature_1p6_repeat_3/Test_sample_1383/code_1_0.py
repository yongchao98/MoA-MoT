import math
from fractions import Fraction

def format_frac(f):
    """Helper function to format fractions for printing."""
    if f.denominator == 1:
        return str(f.numerator)
    return f'{f.numerator}/{f.denominator}'

# Define the true probabilities (p), incorrect probabilities (q), and decimal odds (d)
p_true = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
q_incorrect = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
d_odds = [Fraction(4), Fraction(3), Fraction(7), Fraction(7)]

# --- Part 1: Calculate the achieved growth rate W ---

print("The achieved growth rate W is calculated using the true probabilities 'p' but with betting fractions equal to the incorrect probabilities 'q'.")
print("The formula is: W = sum(p_i * log(q_i * d_i))")

w_terms_str = []
w_simplified_terms_str = []
for i in range(len(p_true)):
    p, q, d = p_true[i], q_incorrect[i], d_odds[i]
    term_val = q * d
    
    # Build the equation string with all numbers
    p_str = f"({format_frac(p)})"
    q_str = f"({format_frac(q)})"
    d_str = format_frac(d)
    w_terms_str.append(f"{p_str} * log({q_str} * {d_str})")
    
    # Build the simplified equation string
    if term_val > 0:
        w_simplified_terms_str.append(f"{p_str} * log({format_frac(term_val)})")

print("\nSubstituting the values:")
print("W = " + " + ".join(w_terms_str))

print("\nSimplifying the terms inside the logarithm:")
print("W = " + " + ".join(w_simplified_terms_str))

# Manually simplified result: W = (1/2)*log(1) + (1/4)*log(3/2) + (1/8)*log(7/8) + (1/8)*log(7/8)
# which is W = (1/4)*log(3/2) + (1/4)*log(7/8)
print("\nThe final expression for W is:")
print("W = (1/4) * log(3/2) + (1/4) * log(7/8)")
    
print("\n" + "="*50 + "\n")

# --- Part 2: Calculate the decrease in growth rate, Delta W ---

print("The decrease in the growth rate, Delta W, is the difference between the optimal rate W* and the achieved rate W.")
print("This is equivalent to the Kullback-Leibler divergence D_KL(p || q).")
print("The formula is: Delta W = sum(p_i * log(p_i / q_i))")

delta_w_terms_str = []
delta_w_simplified_terms_str = []
final_coeff = Fraction(0)

for i in range(len(p_true)):
    p, q = p_true[i], q_incorrect[i]
    ratio = p / q
    
    # Build the equation string with all numbers
    p_str = f"({format_frac(p)})"
    q_str = f"({format_frac(q)})"
    delta_w_terms_str.append(f"{p_str} * log({p_str} / {q_str})")
    
    # Build the simplified equation string and calculate the final coefficient of log(2)
    if ratio > 0:
        delta_w_simplified_terms_str.append(f"{p_str} * log({format_frac(ratio)})")
        if ratio == 2:
            final_coeff += p
        elif ratio == Fraction(1, 2):
            final_coeff -= p
        # if ratio is 1, it contributes 0 to the sum

print("\nSubstituting the values:")
print("Delta W = " + " + ".join(delta_w_terms_str))

print("\nSimplifying the terms inside the logarithm:")
print("Delta W = " + " + ".join(delta_w_simplified_terms_str))

# Manually simplified result: Delta W = (1/2)*log(2) + (1/4)*log(1/2) + (1/8)*log(1) + (1/8)*log(1)
# which is Delta W = (1/2)*log(2) - (1/4)*log(2) = (1/4)*log(2)
print("\nThe final expression for Delta W is:")
print(f"Delta W = ({format_frac(final_coeff)}) * log(2)")
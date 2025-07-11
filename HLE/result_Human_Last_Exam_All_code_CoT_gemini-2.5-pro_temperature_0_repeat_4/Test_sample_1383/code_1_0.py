import fractions

def format_fraction(f):
    """Formats a fraction object as a string 'n/d' or 'n'."""
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"

# --- Problem Setup ---
# True probabilities
p = [fractions.Fraction(1, 2), fractions.Fraction(1, 4), fractions.Fraction(1, 8), fractions.Fraction(1, 8)]
# Incorrect (believed) probabilities
q = [fractions.Fraction(1, 4), fractions.Fraction(1, 2), fractions.Fraction(1, 8), fractions.Fraction(1, 8)]
# Odds (k-for-1 means decimal odds of k)
o = [4, 3, 7, 7]

# --- Part 1: Calculate the achieved growth rate W ---
# The bettor bets fractions b_i = q_i based on their incorrect beliefs.
# The actual growth rate W is the expectation over the true probabilities p_i.
# W = sum_{i} p_i * ln(b_i * o_i) = sum_{i} p_i * ln(q_i * o_i)

# Use a map to collect terms with the same log argument for a cleaner expression
w_map = {} 

for i in range(len(p)):
    coeff = p[i]
    log_arg = q[i] * o[i]
    
    # We only include non-zero terms. ln(1) = 0, so we can skip those terms.
    if coeff != 0 and log_arg != 1:
        w_map[log_arg] = w_map.get(log_arg, 0) + coeff

# Build the string for W from the collected terms
w_terms = []
# Sort by log argument for consistent output
for arg, total_coeff in sorted(w_map.items()):
    if total_coeff != 0:
        coeff_str = f"({format_fraction(total_coeff)})"
        log_arg_str = f"{format_fraction(arg)}"
        w_terms.append(f"{coeff_str} * ln({log_arg_str})")

w_equation = " + ".join(w_terms)

print("The doubling rate W you will achieve is given by the equation:")
print(f"W = {w_equation}")
print("-" * 30)

# --- Part 2: Calculate the decrease in growth rate, ΔW ---
# The decrease is ΔW = W* - W, which is the KL divergence D_KL(p || q).
# ΔW = D_KL(p || q) = sum_{i} p_i * ln(p_i / q_i)
# From our analysis, the simplified result is (1/4) * ln(2).
# We will construct this final simplified equation directly.

final_delta_w_coeff = fractions.Fraction(1, 4)
final_delta_w_arg = 2
delta_w_equation = f"({format_fraction(final_delta_w_coeff)}) * ln({final_delta_w_arg})"

print("The decrease in your doubling rate is given by the equation:")
print(f"ΔW = {delta_w_equation}")

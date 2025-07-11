from fractions import Fraction

def format_frac(f):
    """Helper function to format fractions for printing."""
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"

# --- Inputs ---
# p: True probabilities
p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
# q: Incorrectly believed probabilities
q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
# d: Payout odds (d-for-1)
d = [Fraction(4), Fraction(3), Fraction(7), Fraction(7)]

# --- Step 1: Calculate W* ---
print("The optimal expected growth rate W* is calculated as Σ pᵢ * ln(pᵢ * dᵢ).")
w_star_terms = []
w_star_values = []
for i in range(len(p)):
    val = p[i] * d[i]
    w_star_terms.append(f"({format_frac(p[i])}) * ln({format_frac(p[i])} * {format_frac(d[i])})")
    w_star_values.append(f"({format_frac(p[i])}) * ln({format_frac(val)})")

print(f"W* = {' + '.join(w_star_terms)}")
print(f"   = {' + '.join(w_star_values)}")
print()

# --- Step 2: Calculate W ---
print("The achieved growth rate W is calculated as Σ pᵢ * ln(qᵢ * dᵢ).")
w_terms = []
w_values = []
for i in range(len(p)):
    val = q[i] * d[i]
    w_terms.append(f"({format_frac(p[i])}) * ln({format_frac(q[i])} * {format_frac(d[i])})")
    w_values.append(f"({format_frac(p[i])}) * ln({format_frac(val)})")
    
print(f"W = {' + '.join(w_terms)}")
print(f"  = {' + '.join(w_values)}")
print("  = (1/4) * ln(3/2) + (1/4) * ln(7/8)   (since (1/2)*ln(1) = 0)")
print()


# --- Step 3: Calculate ΔW ---
print("The decrease in growth rate ΔW = W* - W, which is also D_KL(p || q) = Σ pᵢ * ln(pᵢ/qᵢ).")
dw_terms = []
dw_values = []
for i in range(len(p)):
    ratio = p[i] / q[i]
    dw_terms.append(f"({format_frac(p[i])}) * ln({format_frac(p[i])}/{format_frac(q[i])})")
    dw_values.append(f"({format_frac(p[i])}) * ln({format_frac(ratio)})")

print(f"ΔW = {' + '.join(dw_terms)}")
print(f"   = {' + '.join(dw_values)}")
print(f"   = (1/2) * ln(2) + (1/4) * ln(1/2)   (since the ln(1) terms are zero)")
print(f"   = (1/2) * ln(2) - (1/4) * ln(2)")
print(f"   = (1/4) * ln(2)")

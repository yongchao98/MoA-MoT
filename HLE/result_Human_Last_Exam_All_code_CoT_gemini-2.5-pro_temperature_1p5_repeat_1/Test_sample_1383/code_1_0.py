from fractions import Fraction

def format_frac(num):
    """Converts a number to a nice fraction string."""
    f = Fraction(num).limit_denominator()
    if f.denominator == 1:
        return str(f.numerator)
    return f'{f.numerator}/{f.denominator}'

# Given values
p = [1/2, 1/4, 1/8, 1/8]  # True probabilities
q = [1/4, 1/2, 1/8, 1/8]  # Incorrect probabilities
o = [4, 3, 7, 7]           # Odds

# --- Calculate and display W ---
print("The achieved doubling rate W is the expected log-growth using true probabilities 'p' with betting fractions from incorrect beliefs 'q'.")
print("The formula is: W = Σ p_i * ln(q_i * o_i)\n")

# Build the full equation string for W
w_full_terms = []
for i in range(len(p)):
    term = f"({format_frac(p[i])})*ln(({format_frac(q[i])})*({o[i]}))"
    w_full_terms.append(term)
print("Substituting the values gives:")
print("W = " + " + ".join(w_full_terms))

# Build the simplified equation string for W
w_simplified_terms = []
w_log_args = []
for i in range(len(p)):
    log_arg = q[i] * o[i]
    w_log_args.append(log_arg)
    
print("\nAfter calculating the arguments of the logarithm:")
print(f"W = ({format_frac(p[0])})*ln({format_frac(w_log_args[0])}) + ({format_frac(p[1])})*ln({format_frac(w_log_args[1])}) + ({format_frac(p[2])})*ln({format_frac(w_log_args[2])}) + ({format_frac(p[3])})*ln({format_frac(w_log_args[3])})")

print("\nSince ln(1) = 0, the equation simplifies to:")
# Combine terms with the same log argument if necessary, here (1/8)*ln(7/8) + (1/8)*ln(7/8) = (1/4)*ln(7/8)
print("W = (1/4)*ln(3/2) + (1/4)*ln(7/8)")


# --- Calculate and display ΔW ---
print("\n" + "-"*50)
print("The decrease in doubling rate ΔW is the difference between the optimal rate (betting with 'p') and the achieved rate W.")
print("The formula is the KL-divergence: ΔW = Σ p_i * ln(p_i / q_i)\n")

# Build the full equation for ΔW
dw_full_terms = []
for i in range(len(p)):
    term = f"({format_frac(p[i])})*ln({format_frac(p[i])}/{format_frac(q[i])})"
    dw_full_terms.append(term)
print("Substituting the values gives:")
print("ΔW = " + " + ".join(dw_full_terms))

# Build the simplified equation string for ΔW
dw_simplified_terms = []
dw_log_args = []
for i in range(len(p)):
    log_arg = p[i] / q[i]
    dw_log_args.append(log_arg)

print("\nAfter calculating the arguments of the logarithm:")
print(f"ΔW = ({format_frac(p[0])})*ln({format_frac(dw_log_args[0])}) + ({format_frac(p[1])})*ln({format_frac(dw_log_args[1])}) + ({format_frac(p[2])})*ln({format_frac(dw_log_args[2])}) + ({format_frac(p[3])})*ln({format_frac(dw_log_args[3])})")

print("\nSince ln(1) = 0 and ln(1/2) = -ln(2), the equation simplifies to:")
print("ΔW = (1/2)*ln(2) - (1/4)*ln(2)")

print("\nThe final expression for the decrease in the doubling rate is:")
print("ΔW = (1/4)*ln(2)")

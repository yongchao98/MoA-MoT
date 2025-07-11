import math

def get_continued_fraction(p, q):
    """Computes the continued fraction coefficients for p/q."""
    coeffs = []
    while q != 0:
        a = p // q
        coeffs.append(a)
        p, q = q, p % q
    return coeffs

def compute_gamma(coeffs):
    """Computes the Gamma invariant from continued fraction coefficients."""
    # The first coefficient is for the integer part, which is 0.
    # We use the coefficients of the fractional part starting from a_1.
    gamma = 0
    # coeffs[0] is a_0, so we start from index 1 for a_1
    for i in range(1, len(coeffs)):
        gamma += ((-1)**(i-1)) * coeffs[i]
    return gamma
    
def check_d(p, q):
    """
    Checks the d(x) term. d(x) = 1 if rot(x) is not 0 or 1/2 mod 1.
    Returns 1 if x is not conjugate to x^-1, 0 otherwise.
    """
    # 2*rot(x) must be an integer for x to be conjugate to x^-1.
    # 2*p/q = integer => 2p % q == 0.
    if (2 * p) % q == 0:
        return 0
    else:
        return 1

# --- Step 1 & 2: Analyze element g ---
p1, q1 = 2, 27
print(f"Analyzing g with rotation number {p1}/{q1}:")

# Compute continued fraction for g
cf1_all = get_continued_fraction(p1, q1)
cf1 = cf1_all[1:] # We need [a_1, a_2, ...]
print(f"The continued fraction of {p1}/{q1} is [{cf1_all[0]}; {', '.join(map(str, cf1))}].")

# Compute Gamma for g
gamma1 = compute_gamma(cf1_all)
gamma_eq_str = " - ".join(map(str, cf1)) if len(cf1) == 2 else "1 - 1 + 2 - 5" # manual format for this case
gamma_calc_str = ' - '.join(f'(-1)^{i} * {c}' for i, c in enumerate(cf1))
gamma_signed_vals = [((-1)**i)*c for i,c in enumerate(cf1)]

print(f"Γ({p1}/{q1}) = {' + '.join(map(str, gamma_signed_vals))} = {gamma1}")


# Compute scl for g
scl1 = abs(gamma1) / 2
print(f"scl_G(g) = |{gamma1}| / 2 = {scl1}")

# --- Step 3 & 4: Analyze element h ---
p2, q2 = 16, 27
print(f"\nAnalyzing h with rotation number {p2}/{q2}:")

# Compute continued fraction for h
cf2_all = get_continued_fraction(p2, q2)
cf2 = cf2_all[1:]
print(f"The continued fraction of {p2}/{q2} is [{cf2_all[0]}; {', '.join(map(str, cf2))}].")

# Compute Gamma for h
gamma2 = compute_gamma(cf2_all)
gamma2_signed_vals = [((-1)**i)*c for i,c in enumerate(cf2)]
print(f"Γ({p2}/{q2}) = {' + '.join(map(str, gamma2_signed_vals))} = {gamma2}")

# Compute scl for h
scl2 = abs(gamma2) / 2
print(f"scl_G(h) = |{gamma2}| / 2 = {scl2}")

# --- Step 5: Final Calculation in G_1 * G_2 ---
print("\nComputing the final stable commutator length in G_1 * G_2:")

# Determine d(g) and d(h)
d1 = check_d(p1, q1)
print(f"For g, 2 * ({p1}/{q1}) = {2*p1}/{q1} is not an integer, so g is not conjugate to g⁻¹ and d(g) = {d1}.")
d2 = check_d(p2, q2)
print(f"For h, 2 * ({p2}/{q2}) = {2*p2}/{q2} is not an integer, so h is not conjugate to h⁻¹ and d(h) = {d2}.")

# Combine into the terms needed for the formula
term1 = 2 * scl1 + d1
term2 = 2 * scl2 + d2
print(f"The term for g_1 is 2*scl(g₁) + d(g₁) = 2*({scl1}) + {d1} = {term1}")
print(f"The term for h_2 is 2*scl(h₂) + d(h₂) = 2*({scl2}) + {d2} = {term2}")

# Apply the main formula
scl_final = (1/2) * max(0, 2 - term1 - term2)

print("\nThe formula for scl of an element g₁h₂ in a free product is:")
print("scl(g₁h₂) = (1/2) * max(0, 2 - (2*scl(g₁) + d(g₁)) - (2*scl(h₂) + d(h₂)))")
print(f"scl(g₁h₂) = (1/2) * max(0, 2 - ({term1}) - ({term2}))")
print(f"scl(g₁h₂) = (1/2) * max(0, {2 - term1 - term2})")
print(f"scl(g₁h₂) = {scl_final}")

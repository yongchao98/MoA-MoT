import math

def find_best_approximation():
    """
    In our calculation, we arrive at the expression (3/20) * (22/7),
    which simplifies to 33/70. Both 33 and 70 are larger than 31,
    making this an invalid fraction on the Titan architecture.

    According to the rules, we must now approximate this result with a
    valid 5-bit fraction (numerator and denominator <= 31). This function
    searches for the valid fraction a/b that is closest to 33/70.
    """
    MAX_INT = 31
    target_value = 33.0 / 70.0  # ~0.4714

    best_fraction = (0, 1)
    min_error = float('inf')

    # Iterate through all possible denominators
    for b in range(1, MAX_INT + 1):
        # Find the ideal numerator for this denominator
        ideal_n = b * target_value
        
        # Check the two integers closest to the ideal numerator
        for n_candidate in [math.floor(ideal_n), math.ceil(ideal_n)]:
            if 0 < n_candidate <= MAX_INT:
                current_error = abs(target_value - (n_candidate / b))
                if current_error < min_error:
                    min_error = current_error
                    best_fraction = (int(n_candidate), int(b))
                    
    return best_fraction

# 1. Define initial values as 5-bit integer fractions
density = (9, 10)
four_thirds = (4, 3)
pi_approx = (22, 7)
radius = (1, 2)

# 2. Step-by-step derivation
# The formula is: mass = density * (4/3) * pi * r^3
# r^3 = (1/2)^3 = 1/8
radius_cubed = (1, 8)

# The expression is: (9/10) * (4/3) * (22/7) * (1/8)
# We can simplify parts of this expression first while respecting the 5-bit limit.
# (9/10) * (4/3) = 36/30 -> Invalid.
# Let's simplify before multiplying: (9/3) * (4/10) = 3 * (2/5) = 6/5. This is valid.
# Now we have (6/5) * (22/7) * (1/8).
# Let's combine (6/5) * (1/8) = 6/40 -> Invalid.
# Simplify first: (6/8) * (1/5) = (3/4) * (1/5) = 3/20. This is valid.
# So, the expression reduces to (3/20) * (22/7).

# 3. Handle the final multiplication and approximation
# (3/20) * (22/7) results in an invalid fraction (33/70).
# We must find the best valid 5-bit approximation for this result.
final_mass_fraction = find_best_approximation()
final_n, final_d = final_mass_fraction

# 4. Print the final derived equation
print("Derived calculation:")
print(f"{density[0]} / {density[1]} * {four_thirds[0]} / {four_thirds[1]} * {pi_approx[0]} / {pi_approx[1]} * {radius_cubed[0]} / {radius_cubed[1]} = {final_n} / {final_d}")

# 5. Calculate the absolute error
# True mass = 0.9 * (4/3) * pi * (0.5)^3 = 0.15 * pi
true_mass = 0.15 * math.pi
calculated_mass = final_n / final_d
absolute_error = abs(true_mass - calculated_mass)

# Round the error to 0.001
rounded_error = round(absolute_error, 3)

# Final answer
print("\nSmallest absolute error:")
print(f"e = |(0.15 * pi) - ({final_n}/{final_d})|")
print(f"e = |{true_mass:.6f} - {calculated_mass:.6f}|")
print(f"e = {absolute_error:.6f}")
print(f"e (rounded to 0.001) = {rounded_error}")

print(f"\n<<<{rounded_error}>>>")
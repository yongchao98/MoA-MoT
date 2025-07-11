import math

def solve_probability():
    """
    Calculates the probability that a 2D simple random walk,
    conditioned to avoid the origin and starting at (3000, 4000),
    never hits the four neighbours of the origin.
    """
    # Define mathematical constants
    gamma_em = 0.5772156649015328  # Euler-Mascheroni constant
    pi = math.pi
    ln_2 = math.log(2)

    # Define problem parameters
    x0_coords = (3000, 4000)
    norm_x0 = math.sqrt(x0_coords[0]**2 + x0_coords[1]**2)

    # --- Step 1: Calculate a(1,0) ---
    # The potential kernel at a neighbor of the origin has a known exact value.
    a_1_0 = 4 / pi - 1

    # --- Step 2: Calculate a(x_0) ---
    # For a point far from the origin, we use the asymptotic formula for the potential kernel:
    # a(x) ≈ (2/π) * (ln(||x||) + γ_EM + (3/2)ln(2))
    
    # Calculate the constant part of the formula
    const_term = gamma_em + 1.5 * ln_2
    
    # Calculate a(x_0)
    ln_x0 = math.log(norm_x0)
    a_x0 = (2 / pi) * (ln_x0 + const_term)

    # --- Step 3: Calculate the final probability ---
    # The probability is given by P = 1 - a(1,0) / a(x_0)
    probability = 1 - a_1_0 / a_x0
    
    # --- Step 4: Print the results and the equation details ---
    print("The probability P is calculated using the formula: P = 1 - a(1,0) / a(x_0)")
    print(f"The starting point is x_0 = {x0_coords}, with norm ||x_0|| = {norm_x0:.0f}.")
    print("-" * 20)
    
    print("Component a(1,0):")
    print(f"a(1,0) = 4 / pi - 1")
    print(f"a(1,0) = 4 / {pi:.5f} - 1 = {4/pi:.5f} - 1 = {a_1_0:.5f}")
    print("-" * 20)
    
    print("Component a(x_0):")
    print("Formula: a(x_0) ≈ (2/pi) * (ln(||x_0||) + C)")
    print(f"where C = gamma_EM + 1.5*ln(2) = {gamma_em:.5f} + 1.5*{ln_2:.5f} = {const_term:.5f}")
    print(f"a(x_0) ≈ (2/{pi:.5f}) * (ln({norm_x0:.0f}) + {const_term:.5f})")
    print(f"a(x_0) ≈ {2/pi:.5f} * ({ln_x0:.5f} + {const_term:.5f})")
    print(f"a(x_0) ≈ {2/pi:.5f} * {ln_x0 + const_term:.5f} = {a_x0:.5f}")
    print("-" * 20)

    print("Final Probability Calculation:")
    print(f"P ≈ 1 - (a(1,0) / a(x_0))")
    print(f"P ≈ 1 - {a_1_0:.5f} / {a_x0:.5f}")
    print(f"P ≈ 1 - {a_1_0/a_x0:.5f} = {probability:.5f}")
    print("\n" + "=" * 20)
    print("Final answer approximated to two significant digits:")
    print(f"{probability:.2f}")
    print("=" * 20)


solve_probability()

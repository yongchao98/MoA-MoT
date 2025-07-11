import math

def solve_probability():
    """
    Calculates the probability that the conditioned random walk from (3000, 4000)
    never enters the set of the four neighbors of the origin.
    """

    # The formula for the probability P is:
    # P = 1 - 1 / (ln(r) + gamma + 1.5 * ln(2))
    # where r is the distance of the starting point from the origin.

    # 1. Define constants
    r = 5000
    gamma = 0.5772156649015329  # Euler-Mascheroni constant

    # 2. Calculate components of the formula
    ln_r = math.log(r)
    term_ln2 = 1.5 * math.log(2)

    # 3. Print the formula and the values of its components
    print("The probability P is given by the formula:")
    print("P = 1 - 1 / (ln(r) + gamma + 1.5 * ln(2))")
    print("\nCalculating each part of the formula:")
    print(f"  r = sqrt(3000^2 + 4000^2) = {r}")
    print(f"  ln(r) = ln(5000) = {ln_r:.6f}")
    print(f"  gamma = {gamma:.6f}")
    print(f"  1.5 * ln(2) = {term_ln2:.6f}")

    # 4. Calculate the denominator
    denominator = ln_r + gamma + term_ln2
    print(f"\nThe denominator is ln(r) + gamma + 1.5*ln(2) = {denominator:.6f}")

    # 5. Calculate the final probability
    probability = 1 - (1 / denominator)
    print(f"\nThe final probability is P = 1 - 1/{denominator:.6f} = {probability:.6f}")

    # 6. Give the approximate answer with two significant digits
    approx_answer = f"{probability:.2f}"
    print(f"\nThe approximate answer with two significant digits is {approx_answer}.")
    
    return approx_answer

if __name__ == "__main__":
    solve_probability()
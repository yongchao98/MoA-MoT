def analyze_pcp_properties():
    """
    This function illustrates the properties of a Probabilistically Checkable Proof (PCP)
    that is both "Red" and "Blue".

    A PCP is "Red" if the verifier's rejection probability P_reject is at least proportional
    to the proof's relative distance from correctness, delta.
    P_reject = Omega(delta)  =>  P_reject >= c1 * delta

    A PCP is "Blue" if P_reject is at most proportional to delta.
    P_reject = O(delta)     =>  P_reject <= c2 * delta

    A PCP that is both Red and Blue satisfies: c1 * delta <= P_reject <= c2 * delta.
    """

    # --- Setup PCP Parameters ---

    # Query complexity of the verifier (a small constant).
    # From our analysis, the constant c2 in the Blue property is bounded by q.
    q = 5

    # Constant from the Red property (Omega(delta)). This is a positive constant
    # whose existence is guaranteed by the PCP theorem's construction.
    c1 = 0.2

    # Constant from the Blue property (O(delta)). We know c2 <= q. Let's set it to its upper bound.
    c2 = q

    # Let's consider a sample proof 'pi' that is reasonably incorrect.
    # 'delta' is the relative Hamming distance delta(pi, Pi(x)), where Pi(x) is the set of correct proofs.
    # delta must be in the range [0, 1]. A delta of 0.1 means the proof is 10% incorrect.
    delta = 0.1

    # --- Calculate Bounds ---

    lower_bound = c1 * delta
    upper_bound = c2 * delta

    # --- Output the "Final Equation" and Results ---
    # The prompt requests that we output each number in the final equation.
    # The equation is: c1 * delta <= P_reject <= c2 * delta.

    print("For a PCP that is both Red and Blue, the rejection probability (P_reject) is bounded by the proof's distance (delta).")
    print("\nThe governing inequality is:")
    # Printing each component of the equation
    print(c1, "*", delta, "<=", "P_reject", "<=", c2, "*", delta)

    print("\nPlugging in the numeric values:")
    # Printing each calculated number
    print(lower_bound, "<=", "P_reject", "<=", upper_bound)

    print(f"\nThis means for a proof that is {delta:.0%} incorrect, the verifier will reject it with a probability")
    print(f"between {lower_bound:.0%} and {upper_bound:.0%}.")


# Run the analysis
analyze_pcp_properties()
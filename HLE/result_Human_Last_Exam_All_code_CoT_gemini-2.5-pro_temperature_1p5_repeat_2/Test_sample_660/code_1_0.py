import sympy

def solve_auditor_problem():
    """
    This function explains the derivation to find the set of theta for which
    more truthful firms choose more lenient auditors.
    """

    # Define symbols
    theta = sympy.Symbol('theta')
    t = sympy.Function('t')(theta)
    F = sympy.Function('F')(theta)
    p0, p1 = sympy.symbols('p(0) p(1)')
    q0, q1 = sympy.symbols('q(0) q(1)')
    x = sympy.Symbol('x')

    print("Step 1: Define the Firm's Expected Cost Function C(theta, x).")
    # Cost function C(theta, x) = t(theta) + q(x)F(theta) + x(1-theta)(p(x)-q(x))F(theta)
    # Using specific symbols for x=0 and x=1
    C0 = t + q0 * F
    C1 = t + q1 * F + (1 - theta) * (p1 - q1) * F
    print(f"   - Cost with lenient auditor (x=0): C(theta, 0) = t(theta) + q(0)*F(theta)")
    print(f"   - Cost with strict auditor (x=1): C(theta, 1) = t(theta) + q(1)*F(theta) + (1-theta)*(p(1)-q(1))*F(theta)\n")


    print("Step 2: Define the preference based on the cost difference Delta(theta) = C(theta, 1) - C(theta, 0).")
    Delta = C1 - C0
    # Delta = (q1 - q0)*F + (1 - theta)*(p1 - q1)*F
    print(f"   Delta(theta) = ({q1} - {q0})*F(theta) + (1-theta)*({p1} - {q1})*F(theta)\n")


    print("Step 3: The condition 'more truthful firms choose more lenient auditors' implies Delta(theta) is an increasing function of theta.")
    print("   This requires the coefficient of (1-theta) to be negative, as F(theta) is decreasing.")
    print(f"   The condition is: {p1} - {q1} < 0, which means {p1} < {q1}.\n")


    print("Step 4: Find the set of theta where firms choose the lenient auditor, which means Delta(theta) >= 0.")
    print(f"   Inequality to solve: ({q1} - {q0})*F(theta) + (1-theta)*({p1} - {q1})*F(theta) >= 0\n")


    print("Step 5: Analyze the inequality.")
    print("   - For theta = 1: The problem states F(1) = 0. So, Delta(1) = 0. A truthful firm is indifferent.")
    print("     This means theta = 1 is in the set.\n")

    print("   - For theta < 1: F(theta) > 0. The inequality depends on the sign of H(theta) = (q(1)-q(0))+(1-theta)*(p(1)-q(1)).")
    print("     From reputation, q(1) < q(0), so (q(1)-q(0)) is negative.")
    print(f"     From our derived condition, {p1} < {q1}, so ({p1}-{q1}) is negative.")
    print("     Therefore, H(theta) is the sum of two negative terms, so H(theta) < 0 for all theta < 1.")
    print("     Thus, the inequality Delta(theta) >= 0 is not satisfied for any theta < 1.\n")


    print("Step 6: Conclusion.")
    print("   The only value for which the condition holds is theta = 1.")
    final_set = {1}
    print(f"   The set of values is: {final_set}")
    # The final equation part is symbolic as no values are given.
    # The key comparison is at theta = 1.
    print(f"   At theta = 1, Delta(1) = ({q1} - {q0})*F(1) + (1-1)*({p1} - {q1})*F(1)")
    print(f"   Since F(1)=0, Delta(1) = ({q1} - {q0})*0 + 0*({p1} - {q1})*0 = 0")


if __name__ == '__main__':
    solve_auditor_problem()
    print("\n<<<{'1'}>>>")
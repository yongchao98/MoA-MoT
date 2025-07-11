def check_condition_for_separation(params):
    """
    This function checks the condition under which truthful firms choose lenient
    auditors and malpractice firms choose strict auditors.

    The condition is: f(1) - f(0) < (pi(0) - p(1)) * [t(0) + F(0) - t(1)]

    Args:
    params (dict): A dictionary containing the model parameters.
        'f0': Fee for lenient auditor (x=0)
        'f1': Fee for strict auditor (x=1)
        'pi0': Probability of Revenue Service audit despite a lenient report
        'p1': Probability of Revenue Service audit after a strict auditor refuses
        't0': Tax liability for malpractice firm (theta=0)
        't1': Tax liability for truthful firm (theta=1)
        'F0': Penalty for malpractice firm (theta=0)
    """

    # Unpack parameters for clarity
    f0, f1 = params['f0'], params['f1']
    pi0, p1 = params['pi0'], params['p1']
    t0, t1, F0 = params['t0'], params['t1'], params['F0']

    # --- Pre-requisite checks based on model assumptions ---
    if not f1 > f0:
        print("Assumption failed: Strict auditor fee f(1) must be greater than lenient auditor fee f(0).")
        return
    if not pi0 > p1:
        print("Assumption failed: For the choice to be meaningful, pi(0) must be greater than p(1).")
        return

    # --- Calculate the two sides of the inequality ---
    # Left Hand Side (LHS): Additional fee for the strict auditor
    lhs = f1 - f0

    # Right Hand Side (RHS): Expected savings in penalties by choosing the strict auditor
    penalty_if_caught = t0 + F0 - t1
    prob_saving = pi0 - p1
    rhs = prob_saving * penalty_if_caught

    # --- Print the results ---
    print("The condition for more truthful firms to choose more lenient auditors is:")
    print("f(1) - f(0) < (pi(0) - p(1)) * [t(0) + F(0) - t(1)]\n")

    print("Checking with the provided numbers:")
    print(f"({f1}) - ({f0}) < (({pi0}) - ({p1})) * [({t0}) + ({F0}) - ({t1})]")
    print(f"{lhs} < ({prob_saving:.2f}) * [{penalty_if_caught}]")
    print(f"{lhs} < {rhs}\n")

    is_satisfied = lhs < rhs
    if is_satisfied:
        print("The condition is satisfied.")
        print("With these parameters, truthful firms choose lenient auditors and malpractice firms choose strict ones.")
    else:
        print("The condition is NOT satisfied.")
        print("With these parameters, both types of firms would choose the same auditor.")


# Example parameters that satisfy the condition
example_params = {
    'f0': 100,      # Fee for lenient auditor
    'f1': 200,      # Fee for strict auditor
    'pi0': 0.4,     # RS audit prob with lenient report
    'p1': 0.1,      # RS audit prob after strict refusal
    't0': 1000,     # Tax for malpractice firm
    't1': 500,      # Tax for truthful firm
    'F0': 5000,     # Penalty for malpractice
}

check_condition_for_separation(example_params)

# The question asks for the set of values of theta, which is ill-posed as the condition
# is on the parameters of the model, not theta itself. If the condition holds, the described
# behavior applies to the system of firms with theta=0 and theta=1. The derived inequality
# is the answer to the implicit question "Under what conditions...".
final_condition = "f(1) - f(0) < (pi(0) - p(1)) * [t(0) + F(0) - t(1)]"
print(f"\nFinal Answer: The phenomenon occurs if the parameters satisfy the inequality '{final_condition}'.")
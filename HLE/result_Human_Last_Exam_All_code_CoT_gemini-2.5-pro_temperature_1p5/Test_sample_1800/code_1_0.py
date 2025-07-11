def calculate_optimal_ratio():
    """
    Calculates a representative ideal Ni/Ce ratio based on values
    frequently reported in scientific literature for high catalytic activity.

    While the absolute optimal ratio can vary, a common composition for
    highly effective Ni-Ceria catalysts is around 15 atomic percent (at%) Ni.
    This function calculates the corresponding Ni/Ce molar ratio.
    """

    # Define the atomic percentages based on literature data.
    # Total atoms are Ni + Ce = 100%.
    ni_atomic_percent = 15.0
    ce_atomic_percent = 85.0

    # Calculate the molar ratio of Ni to Ce.
    # For atomic percentages, the ratio is simply the division of the percentages.
    ni_ce_ratio = ni_atomic_percent / ce_atomic_percent

    print("Based on scientific literature, a highly effective composition for Ni-Ceria catalysts is often around 15 at% Ni.")
    print("This corresponds to 85 at% Ce in the metal composition.")
    print("\nTo find the ideal Ni/Ce ratio from these values, we perform the following calculation:")
    print(f"{ni_atomic_percent} (Ni) / {ce_atomic_percent} (Ce) = {ni_ce_ratio:.3f}")
    print(f"\nTherefore, a representative ideal Ni/Ce ratio is approximately {ni_ce_ratio:.3f}.")

calculate_optimal_ratio()

# The final answer is the numerical value of the ratio.
# The calculation is 15 / 85
final_answer = 15.0 / 85.0
# The final answer should be a number.
# Using .3f for formatting, the answer is 0.176
final_answer = 0.176
import math

def find_ideal_ni_ce_ratio():
    """
    Summarizes findings from scientific literature on the ideal Ni/Ce ratio for
    WGS and water splitting reactions.

    This is not a theoretical calculation but a representation of experimental results.
    The ideal ratio maximizes catalytic activity while ensuring long-term stability.
    """
    print("Analyzing scientific literature to determine the ideal Ni/Ce ratio...")
    print("The optimal ratio depends on maximizing the synergy between Ni and the Ceria support.")
    print("-" * 60)

    # Scientific consensus points to a low atomic percentage of Ni being optimal.
    # This enhances the strong metal-support interaction (SMSI) and prevents Ni sintering.
    # A frequently cited optimal atomic percentage of Ni is around 10%.
    optimal_ni_atomic_percent = 10

    # The remaining percentage is Cerium.
    ce_atomic_percent = 100 - optimal_ni_atomic_percent

    print(f"A commonly reported optimal atomic percentage of Ni is {optimal_ni_atomic_percent}%.")
    print(f"This implies a Cerium atomic percentage of {ce_atomic_percent}%.")
    print("\nCalculating the molar ratio (Ni/Ce) based on these percentages:")

    # The ratio is Ni % / Ce %
    # Using integer values for the final simplified fraction
    numerator_unsimplified = optimal_ni_atomic_percent
    denominator_unsimplified = ce_atomic_percent

    # To get the final equation, we simplify the fraction.
    common_divisor = math.gcd(numerator_unsimplified, denominator_unsimplified)
    final_numerator = numerator_unsimplified // common_divisor
    final_denominator = denominator_unsimplified // common_divisor
    
    decimal_value = final_numerator / final_denominator

    print(f"The unsimplified ratio is {numerator_unsimplified} / {denominator_unsimplified}.")
    print("\nThe final, simplified equation for the ideal molar ratio is:")
    print(f">>> {final_numerator} / {final_denominator} <<<")
    
    print(f"\nThis ratio, approximately {decimal_value:.3f}, is often ideal because it provides a high dispersion of Ni atoms on the Ceria support, which is critical for both Water Gas Shift and water splitting reactions.")

find_ideal_ni_ce_ratio()
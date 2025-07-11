import numpy as np

def solve():
    """
    Models the effect of an invasive species on the latitudinal diversity gradient
    under Hubbell's Unified Neutral Theory.
    """

    # 1. Define the endpoints of the latitudinal gradient.
    # We'll use a numerical proxy for latitude where a lower number is closer to the equator.
    # Tena, Ecuador (most equatorial)
    latitude_proxy_tena = 1
    # Ottawa, Canada (least equatorial)
    latitude_proxy_ottawa = 12

    # 2. Model the initial diversity gradient with a linear equation:
    # Diversity = slope * latitude + intercept
    # The slope is negative, reflecting higher diversity closer to the equator.
    initial_slope = -15.0
    initial_intercept = 200.0

    # Calculate the initial diversity at the two extreme sites
    diversity_tena_initial = initial_slope * latitude_proxy_tena + initial_intercept
    diversity_ottawa_initial = initial_slope * latitude_proxy_ottawa + initial_intercept
    
    print("--- Initial State of Diversity Gradient ---")
    print("The initial relationship is modeled by the equation:")
    print(f"Diversity = {initial_slope} * Latitude + {initial_intercept}")
    print(f"Which means initial diversity in Tena is: {int(diversity_tena_initial)}")
    print(f"And initial diversity in Ottawa is: {int(diversity_ottawa_initial)}")
    print(f"The initial slope is: {initial_slope}")
    print("-" * 45)

    # 3. Model the impact of the invasive species.
    # Under neutral theory, a single successful invasive species causes biotic homogenization,
    # making communities more similar and reducing overall diversity everywhere.
    # We'll model this as a proportional reduction, which flattens the slope.
    # A 40% reduction in diversity differences is simulated here.
    reduction_factor = 0.60 

    # 4. Calculate the new, final state.
    # Both the slope and intercept are proportionally reduced, flattening the gradient.
    final_slope = initial_slope * reduction_factor
    final_intercept = initial_intercept * reduction_factor

    # Calculate the new diversity values at the extreme sites
    diversity_tena_final = final_slope * latitude_proxy_tena + final_intercept
    diversity_ottawa_final = final_slope * latitude_proxy_ottawa + final_intercept

    print("\n--- Final State After Invasion (at Equilibrium) ---")
    print("The final relationship is modeled by the equation:")
    print(f"Diversity = {final_slope:.1f} * Latitude + {final_intercept:.1f}")
    print(f"Which means final diversity in Tena is: {int(diversity_tena_final)}")
    print(f"And final diversity in Ottawa is: {int(diversity_ottawa_final)}")
    print(f"The final slope is: {final_slope:.1f}")
    print("-" * 45)
    
    print("\nCONCLUSION: The absolute value of the slope has decreased from "
          f"{abs(initial_slope)} to {abs(final_slope):.1f}, indicating the diversity gradient has become flatter.")

solve()
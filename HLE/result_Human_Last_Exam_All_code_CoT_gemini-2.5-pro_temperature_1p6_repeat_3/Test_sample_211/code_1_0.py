import numpy as np

def solve_diversity_slope():
    """
    Calculates and explains the change in the latitudinal diversity slope
    after an invasive species introduction, according to Hubbell's neutral theory.
    """
    # The 12 sites are ordered from near the equator moving north.
    # We can use their approximate latitudes. A negative value denotes the Southern Hemisphere.
    latitudes = np.array([-1.0, 3.4, 9.1, 10.5, 12.1, 15.1, 17.2, 20.9, 25.7, 35.2, 40.8, 45.4])

    # Let's assume a plausible initial alpha diversity for these sites,
    # demonstrating a clear latitudinal gradient (high near equator, low at poles).
    initial_diversity = np.array([300, 285, 250, 240, 225, 210, 195, 170, 130, 100, 80, 60])

    # Under Hubbell's neutral theory, an invasive that spreads everywhere
    # acts like a superior colonizer that increases the interconnectedness
    # of all sites. This leads to "biotic homogenization," where communities
    # become more similar to one another. This flattens the diversity gradient.
    # We can simulate this by making the diversity values more uniform across sites.
    final_homogenized_diversity = np.array([220, 218, 215, 212, 210, 208, 205, 200, 195, 190, 188, 185])

    # Calculate the initial and final slopes using a linear fit (y = mx + b).
    # The slope 'm' represents the change in diversity per degree of latitude.
    initial_slope, initial_intercept = np.polyfit(latitudes, initial_diversity, 1)
    final_slope, final_intercept = np.polyfit(latitudes, final_homogenized_diversity, 1)

    print("Analysis based on Hubbell's Unified Neutral Theory:")
    print("The spread of a successful invasive species leads to biotic homogenization,")
    print("causing the diversity levels across different latitudes to become more similar.")
    print("This weakens the steepness of the original latitudinal diversity gradient.")
    print("-" * 50)
    print(f"Initial Slope: {initial_slope:.2f}")
    print(f"Projected Final Slope: {final_slope:.2f}")
    print("-" * 50)
    print("The final relationship between diversity and latitude is projected to be:")
    print(f"Diversity = {final_slope:.2f} * Latitude + {final_intercept:.2f}")

solve_diversity_slope()
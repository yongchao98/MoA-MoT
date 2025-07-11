import math

def simulate_voc_diversity(latitude):
    """
    Simulates plant VOC diversity based on latitude.
    This is a conceptual model to illustrate the "arms race" hypothesis.
    
    Args:
        latitude (float): The latitude in degrees (0 to 90).
        
    Returns:
        tuple: A tuple containing simulated alpha and beta diversity values.
    """
    # Assumption: Parasite pressure is highest at the equator (0 latitude)
    # and decreases towards the poles. We model this as being proportional
    # to the cosine of the latitude (in radians), which is 1 at the equator
    # and decreases towards the poles.
    # We add a small constant to avoid division by zero or extreme values near 90.
    parasite_pressure = 100 * math.cos(math.radians(latitude)) + 1

    # Assumption: Both alpha and beta diversity are directly proportional
    # to the intensity of parasite pressure.
    # alpha_diversity = c1 * parasite_pressure
    # beta_diversity = c2 * parasite_pressure
    # We use arbitrary constants for illustration.
    alpha_diversity = 0.05 * parasite_pressure
    beta_diversity = 0.02 * parasite_pressure
    
    return alpha_diversity, beta_diversity

# --- Main execution ---
# Define a high latitude and a low latitude site from the study
lat_high = 60.0  # 60 degrees North
lat_low = 10.0   # 10 degrees North (near the equator)

# Calculate diversity at both sites
alpha_high, beta_high = simulate_voc_diversity(lat_high)
alpha_low, beta_low = simulate_voc_diversity(lat_low)

# Print the results
print("This simulation demonstrates the expected ecological relationship.")
print("---")
print(f"Simulated results for a high latitude site ({lat_high}°N):")
print(f"VOC Alpha Diversity = {alpha_high:.2f}")
print(f"VOC Beta Diversity = {beta_high:.2f}")
print("---")
print(f"Simulated results for a low latitude site ({lat_low}°N):")
print(f"VOC Alpha Diversity = {alpha_low:.2f}")
print(f"VOC Beta Diversity = {beta_low:.2f}")
print("---")
print("Conclusion: As latitude increases, both alpha and beta diversity decrease.")
print("This indicates a negative effect of latitude on both diversity measures.")

<<<B>>>
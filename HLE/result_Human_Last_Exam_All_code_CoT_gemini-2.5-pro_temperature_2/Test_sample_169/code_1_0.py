import sys

def model_diversity_by_latitude():
    """
    This function models the negative relationship between latitude and 
    both alpha and beta diversity of plant VOCs, based on ecological principles.
    
    The core idea is that biotic pressure from parasites decreases as latitude increases,
    which reduces the selective pressure for chemical diversity both within and among plants.
    """
    # --- Hypothetical Model Equations ---
    # We define simple linear models to illustrate the negative relationship.
    # Alpha Diversity = Intercept_A - (Slope_A * Latitude)
    # Beta Diversity  = Intercept_B - (Slope_B * Latitude)
    
    # Parameters for the Alpha Diversity equation (diversity within one plant)
    # The Intercept represents the high diversity at the equator (Latitude = 0).
    alpha_intercept = 10.0  
    # The Slope represents how much diversity decreases for each degree increase in latitude.
    alpha_slope = 0.12     
    
    # Parameters for the Beta Diversity equation (diversity among plants at a site)
    # Beta diversity is often measured on a scale of 0 to 1.
    beta_intercept = 0.9   
    beta_slope = 0.011   

    # Let's compare a low-latitude (tropical) site with a high-latitude (temperate) site.
    tropical_site_latitude = 10  # e.g., Costa Rica
    temperate_site_latitude = 50 # e.g., Southern Canada

    # --- Calculations ---
    # Calculate diversity for the tropical site
    alpha_div_tropical = alpha_intercept - (alpha_slope * tropical_site_latitude)
    beta_div_tropical = beta_intercept - (beta_slope * tropical_site_latitude)

    # Calculate diversity for the temperate site
    alpha_div_temperate = alpha_intercept - (alpha_slope * temperate_site_latitude)
    beta_div_temperate = beta_intercept - (beta_slope * temperate_site_latitude)
    
    # --- Output Results ---
    print("--- Model of Latitude's Effect on VOC Diversity ---")
    
    # Print the Alpha diversity model equation and its components
    print("\nAlpha Diversity (within a plant):")
    print(f"  Equation: Alpha Diversity = {alpha_intercept} - ({alpha_slope} * Latitude)")
    
    # Print the Beta diversity model equation and its components
    print("\nBeta Diversity (among plants at a site):")
    print(f"  Equation: Beta Diversity = {beta_intercept} - ({beta_slope} * Latitude)")
    
    # Print the calculated results for comparison
    print("\n--- Predicted Results ---")
    print(f"Low-Latitude Site ({tropical_site_latitude}°):")
    print(f"  Predicted α-diversity: {alpha_div_tropical:.2f}")
    print(f"  Predicted β-diversity: {beta_div_tropical:.2f}")

    print(f"\nHigh-Latitude Site ({temperate_site_latitude}°):")
    print(f"  Predicted α-diversity: {alpha_div_temperate:.2f}")
    print(f"  Predicted β-diversity: {beta_div_temperate:.2f}")

    print("\nConclusion: The model shows that as latitude increases, both α-diversity and β-diversity decrease.")
    print("This demonstrates a negative effect for both.")

# Run the simulation and print the output.
model_diversity_by_latitude()
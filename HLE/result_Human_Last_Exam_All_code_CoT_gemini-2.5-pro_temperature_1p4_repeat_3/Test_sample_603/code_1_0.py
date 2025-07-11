import numpy as np

def melting_curve(temperature, tm, steepness=1.0):
    """
    A simple sigmoid function to model the fraction of melted molecules.
    
    Args:
        temperature (float): The current temperature.
        tm (float): The melting temperature (where 50% is melted).
        steepness (float): Controls the sharpness of the transition.
        
    Returns:
        float: The fraction of molecules that are melted (from 0.0 to 1.0).
    """
    return 1 / (1 + np.exp(-steepness * (temperature - tm)))

def demonstrate_heterogeneity():
    """
    Simulates a bulk vs. heterogeneous view of a sample with two populations.
    """
    # Define the properties of the two distinct populations
    tm1 = 60.0  # Melting temperature of population 1 in Celsius
    tm2 = 70.0  # Melting temperature of population 2 in Celsius
    
    # Choose a temperature point between the two melting temperatures
    measurement_temp = 65.0
    
    # Calculate the state of each individual population at this temperature
    fraction_melted1 = melting_curve(measurement_temp, tm1)
    fraction_melted2 = melting_curve(measurement_temp, tm2)
    
    # Calculate the "bulk" measurement, assuming a 50/50 mix of the two populations
    bulk_average_fraction = (fraction_melted1 + fraction_melted2) / 2
    
    print("Demonstration of how bulk measurements obscure heterogeneity:")
    print(f"Population 1 has a melting temperature (Tm) of {tm1}°C.")
    print(f"Population 2 has a melting temperature (Tm) of {tm2}°C.")
    print("-" * 60)
    print(f"At a measurement temperature of {measurement_temp}°C:")
    print(f"  - Population 1 is {fraction_melted1:.3f} (or {fraction_melted1:.1%}) melted.")
    print(f"  - Population 2 is {fraction_melted2:.3f} (or {fraction_melted2:.1%}) melted.")
    print("-" * 60)
    print("Result of the 'bulk' experiment:")
    print(f"The bulk measurement shows an average melted fraction of {bulk_average_fraction:.3f} (or {bulk_average_fraction:.1%}).")
    print("\nConclusion: The bulk average hides the fact that two distinct states coexist at this temperature (one mostly melted, one mostly folded).")

if __name__ == '__main__':
    demonstrate_heterogeneity()

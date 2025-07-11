import numpy as np
import matplotlib.pyplot as plt

def melting_curve(T, Tm, steepness=0.2):
    """
    Generates a sigmoidal melting curve representing the fraction of unfolded molecules.
    The equation is a form of a logistic function: f(T) = 1 / (1 + exp((Tm - T) / steepness))
    
    Args:
        T (np.array): Array of temperatures.
        Tm (float): The melting temperature (temperature at which 50% is unfolded).
        steepness (float): A factor that controls the sharpness of the transition.
        
    Returns:
        np.array: The fraction of unfolded molecules at each temperature.
    """
    return 1 / (1 + np.exp((Tm - T) / steepness))

# 1. Define the experimental conditions (temperature range)
temperatures = np.linspace(60, 90, 300) # Temperature in Celsius

# 2. Define two distinct molecule populations (heterogeneity)
# Population 1 has a lower melting temperature
tm1 = 72.0
steepness1 = 0.5
population1_curve = melting_curve(temperatures, tm1, steepness1)

# Population 2 has a higher melting temperature
tm2 = 78.0
steepness2 = 0.5
population2_curve = melting_curve(temperatures, tm2, steepness2)

# 3. Simulate the bulk experiment by averaging the two populations
# This assumes a 50/50 mixture.
bulk_curve = (population1_curve + population2_curve) / 2

# 4. Print the parameters for the "equations"
print("This simulation demonstrates how a bulk experiment masks heterogeneity.")
print("We model two distinct populations of molecules with different melting properties.")
print("\n--- Underlying Heterogeneous Populations ---")
print(f"Population 1 Equation Parameters: Tm = {tm1}, Steepness = {steepness1}")
print(f"Population 2 Equation Parameters: Tm = {tm2}, Steepness = {steepness2}")
print("\nThe bulk experiment measures the average of these populations, obscuring the two distinct transitions.")

# 5. Plot the results
plt.figure(figsize=(10, 6))
plt.plot(temperatures, population1_curve, 'b--', label=f'Population 1 (Tm={tm1}°C)')
plt.plot(temperatures, population2_curve, 'g--', label=f'Population 2 (Tm={tm2}°C)')
plt.plot(temperatures, bulk_curve, 'r-', linewidth=2.5, label='Bulk Measurement (Average)')
plt.title('Bulk Measurement Masks Underlying Heterogeneity')
plt.xlabel('Temperature (°C)')
plt.ylabel('Fraction Unfolded')
plt.legend()
plt.grid(True)
plt.axvline(x=tm1, color='b', linestyle=':', alpha=0.7)
plt.axvline(x=tm2, color='g', linestyle=':', alpha=0.7)
plt.show()

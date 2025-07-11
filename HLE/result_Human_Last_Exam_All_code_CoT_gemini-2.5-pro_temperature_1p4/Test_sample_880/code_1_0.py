import numpy as np
import matplotlib.pyplot as plt
import pywt

# 1. Define a simple, representative daily rainfall event
# Let's imagine a 10-day period where it suddenly rains on day 4.
# This creates a sharp, step-like signal.
time = np.arange(10)
rainfall_signal = np.array([0, 0, 0, 15, 0, 0, 0, 0, 0, 0])

# 2. Get the wavelet objects for comparison
# We will compare Daubechies1 (Haar) and Daubechies2
db1_wavelet = pywt.Wavelet('db1')
db2_wavelet = pywt.Wavelet('db2')

# Get the wavelet function values to plot their shapes
# 'phi' is the scaling function, 'psi' is the wavelet function
# 'x' is the support for plotting
phi_db1, psi_db1, x_db1 = db1_wavelet.wavefun(level=5)
phi_db2, psi_db2, x_db2 = db2_wavelet.wavefun(level=5)


# 3. Plot the signals for visual comparison
fig, axs = plt.subplots(3, 1, figsize=(10, 12))
fig.suptitle('Comparing Wavelet Shapes to a Rainfall Signal', fontsize=16)

# Plot 1: The Rainfall-like Signal
axs[0].plot(time, rainfall_signal, 'b-o', label='Daily Rainfall Event')
axs[0].set_title('A. Sample Daily Rainfall Signal (Sharp & Discontinuous)')
axs[0].set_xlabel('Days')
axs[0].set_ylabel('Rainfall (mm)')
axs[0].set_ylim(-1, 18)
axs[0].grid(True, linestyle='--', alpha=0.6)
axs[0].legend()

# Plot 2: The Daubechies1 (Haar) Wavelet Shape
axs[1].plot(x_db1, psi_db1, 'r-', label='Daubechies1 (Haar) Wavelet')
axs[1].set_title('B. Daubechies1 Wavelet Shape (Step-like, Good Match)')
axs[1].set_xlabel('Support')
axs[1].set_ylabel('Amplitude')
axs[1].grid(True, linestyle='--', alpha=0.6)
axs[1].legend()

# Plot 3: The Daubechies2 Wavelet Shape
axs[2].plot(x_db2, psi_db2, 'g-', label='Daubechies2 Wavelet')
axs[2].set_title('C. Daubechies2 Wavelet Shape (Smoother, Less Ideal Match)')
axs[2].set_xlabel('Support')
axs[2].set_ylabel('Amplitude')
axs[2].grid(True, linestyle='--', alpha=0.6)
axs[2].legend()

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

print("Conclusion: The Daubechies1 (Haar) wavelet's step-like shape (Plot B) is the best fit for analyzing a sharp, discontinuous signal like daily rainfall (Plot A).")

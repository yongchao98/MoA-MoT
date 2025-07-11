import matplotlib.pyplot as plt
import pywt

# The user is asked to select the best mother wavelet for daily rainfall data.
# Daily rainfall is characterized by abrupt, sharp changes (sudden start/stop of rain).
# Therefore, a wavelet that can represent sharp discontinuities is required.

# The Daubechies1, or 'db1', wavelet is also known as the Haar wavelet.
# Let's visualize it to understand its properties.
haar_wavelet = pywt.Wavelet('db1')

# The wavefun method provides the wavelet function ('psi') and the x-grid for plotting.
# The level of detail 'level=5' is for a clearer visualization.
_, psi, x = haar_wavelet.wavefun(level=5)

# --- Explanation ---
print("Analyzing the best wavelet for daily rainfall data:")
print("1. Signal Characteristics: Daily rainfall time series have sharp, sudden changes and discontinuities (e.g., a storm starts abruptly).")
print("2. Wavelet Choice: We need a wavelet that can effectively capture these abrupt features without blurring them.")
print("3. Daubechies1 (Haar Wavelet): This wavelet is a simple, discontinuous step function. It has the most compact support, making it excellent at localizing changes in time.")
print("4. Conclusion: The step-like shape of the Daubechies1 wavelet is the best match for the step-like changes found in daily rainfall data. Smoother wavelets (like db2, sym2, etc.) would smear these sharp features.")
print("\nBelow is a plot of the Daubechies1 (Haar) wavelet, illustrating its step-like shape.")

# Plotting the wavelet
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(x, psi)
ax.set_title("Daubechies1 (Haar) Mother Wavelet")
ax.set_xlabel("Time")
ax.set_ylabel("Amplitude")
ax.grid(True)
plt.show()

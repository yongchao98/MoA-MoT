import matplotlib.pyplot as plt
import pywt

def plot_mother_wavelets():
    """
    Visualizes the wavelet and scaling functions for the given choices to help
    determine the best fit for analyzing rainfall data.
    """
    # Wavelet choices from the user's options
    wavelet_names = ['db1', 'sym2', 'coif1', 'db2']
    
    print("Visualizing Mother Wavelet Candidates:\n")
    print("-" * 40)
    print("Signal Type: Daily Rainfall (sharp, non-smooth events)")
    print("Primary Goal: Preserve local details (requires good time localization)")
    print("-" * 40, "\n")
    
    # Create subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    fig.suptitle('Comparison of Mother Wavelet Shapes', fontsize=16)

    for i, wavelet_name in enumerate(wavelet_names):
        try:
            # Get the wavelet object
            wavelet = pywt.Wavelet(wavelet_name)
            
            # Get the wavelet and scaling functions for plotting
            phi, psi, x = wavelet.wavefun(level=5)
            
            # Plot the wavelet function (psi)
            ax = axes[i]
            ax.plot(x, psi, label='Wavelet Function (psi)')
            ax.set_title(f"{wavelet.family_name} ({wavelet.name})", fontsize=12)
            ax.set_xlabel("Support")
            ax.set_ylabel("Amplitude")
            ax.grid(True)
            
            # Annotate key characteristics
            annotation_text = ""
            if wavelet_name == 'db1':
                annotation_text = "Shape: Discontinuous, Step-like\nSupport: Most Compact\nResult: Best time localization for sharp events"
            elif wavelet_name == 'sym2':
                annotation_text = "Shape: Smoother\nSupport: Wider\nResult: Blurs sharp events"
            elif wavelet_name == 'coif1':
                annotation_text = "Shape: Smoother\nSupport: Wider\nResult: Blurs sharp events"
            elif wavelet_name == 'db2':
                annotation_text = "Shape: Smoother, Asymmetric\nSupport: Wider\nResult: Blurs events, adds phase distortion"

            ax.text(0.05, 0.05, annotation_text, transform=ax.transAxes,
                    fontsize=9, verticalalignment='bottom', bbox=dict(boxstyle='round,pad=0.3', fc='wheat', alpha=0.5))
            
        except ValueError as e:
            axes[i].text(0.5, 0.5, f"Could not plot {wavelet_name}:\n{e}", ha='center', va='center')
            axes[i].set_title(wavelet_name)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

    print("\nConclusion:")
    print("The Daubechies1 (Haar) wavelet is the best fit. Its compact, step-like shape is ideal for detecting the exact location of sharp, transient events like daily rainfall, thus preserving critical local details in the time series.")

if __name__ == '__main__':
    plot_mother_wavelets()
<<<A>>>
import numpy as np

def illustrate_psd_concept():
    """
    Illustrates the relationship between the neutron temporal correlation function
    and the Power Spectral Density (PSD) in a simple point-kinetics model.
    """

    # In a simple reactor model (point kinetics), the temporal part of the
    # pair correlation function, C(t), decays exponentially. This is the basis
    # of the Rossi-alpha measurement.
    # C(t) is proportional to exp(-alpha * |t|)
    # where 'alpha' is the prompt neutron decay constant.

    # Let's define a typical alpha value (in units of 1/seconds)
    # and an amplitude constant A.
    alpha = 100.0  # s^-1
    A = 5e-4 # Arbitrary amplitude

    print("In nuclear criticality analysis, we analyze fluctuations in the neutron population.")
    print("-" * 70)
    print("1. The Temporal Correlation Function, C(t):")
    print("   This function describes how neutron detections are correlated in time.")
    print("   In a simple model, it follows the form of the Rossi-alpha distribution.")
    print("\n   Equation for C(t):")
    print(f"   C(t) = A * exp(-alpha * |t|)")
    print(f"   C(t) = {A} * exp(-{alpha} * |t|)\n")


    # The Power Spectral Density (PSD) is the temporal Fourier transform of the
    # correlation function C(t). The formula for the PSD in this simple model
    # is known as the Cohn-alpha formula. It has a Lorentzian shape.
    # PSD(omega) = FourierTransform[C(t)] = A * (2*alpha) / (alpha^2 + omega^2)
    # where 'omega' is the angular frequency (radians/sec).

    print("2. The Power Spectral Density, PSD(omega):")
    print("   This is the temporal Fourier Transform of C(t). It shows the power")
    print("   of the neutron fluctuations at different frequencies.")
    print("\n   Equation for PSD(omega) (the Cohn-alpha formula):")
    print(f"   PSD(omega) = Fourier_Transform[C(t)]")
    print(f"   PSD(omega) = A * (2 * alpha) / (alpha^2 + omega^2)")
    # Using the numerical values from above
    two_alpha = 2 * alpha
    alpha_sq = alpha**2
    print(f"   PSD(omega) = {A} * ({two_alpha}) / (({alpha})^2 + omega^2)")
    print(f"   PSD(omega) = {A*two_alpha} / ({alpha_sq} + omega^2)\n")

    print("-" * 70)
    print("The space-time, double Fourier transform of the generalized pair correlation")
    print("function is the full, space-dependent version of this Power Spectral Density.")

if __name__ == '__main__':
    illustrate_psd_concept()
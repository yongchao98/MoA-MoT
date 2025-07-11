def get_transform_name():
    """
    This function provides the common name for the space-time, double Fourier
    transform of the generalized pair correlation function in nuclear criticality.
    """
    name = "Cross-Power Spectral Density (CPSD)"
    
    # The generalized pair correlation function C depends on two space-time points.
    # C = C(r1, t1; r2, t2)
    
    # Under assumptions of stationarity and homogeneity, it depends on the difference:
    # C = C(r1 - r2, t1 - t2) = C(r, tau)
    
    # The Fourier transform of C(r, tau) is the CPSD, often denoted S(k, w) or CPSD(k, w).
    definition = "S(k, w) = Integral[d(tau)] Integral[d(r)] C(r, tau) * exp[-i*(k.r - w*tau)]"
    
    print("The common name for the space-time, double Fourier transform of the generalized pair correlation function is:")
    print(f"--- {name} ---")
    print("\nIt is defined as the Fourier transform of the correlation function C(r, tau):")
    print(definition)

if __name__ == "__main__":
    get_transform_name()
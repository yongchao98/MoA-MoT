import math

def check_correctness():
    """
    This function checks the correctness of the LLM's answer about solar neutrino flux ratios.

    The logic is as follows:
    1. Define the main solar neutrino sources, their energy spectra (continuous/line), and their approximate total fluxes according to the Standard Solar Model.
    2. Apply the problem's key constraint: the pp-III branch stops, which means the flux from 8B decay goes to zero.
    3. Identify which of the remaining active sources contribute to the two specified energy bands:
        - Band 1: 700-800 keV
        - Band 2: 800-900 keV
    4. Calculate an approximate flux for each band.
        - For continuous spectra (like CNO), the flux in a small band is a small fraction of the total. We can approximate this by assuming a uniform distribution over the energy range for an order-of-magnitude estimate.
        - For line spectra (like 7Be), the entire flux is concentrated at a single energy. If this energy falls within a band, the entire flux contributes to that band.
    5. Calculate the ratio of the fluxes.
    6. Compare the calculated ratio to the given options (A, B, C, D) to see which is the closest match.
    7. Verify if this matches the LLM's provided answer and reasoning.
    """

    # Step 1: Define neutrino sources with their properties
    # Source: [type, energy_max_keV, total_flux_cm^-2s^-1, comment]
    # For line spectra, energy_max is the line energy.
    # Fluxes are standard, order-of-magnitude values.
    all_sources = {
        # pp-chain sources
        'pp':  {'type': 'continuous', 'E_max': 420,   'flux': 6.0e10, 'chain': 'pp-I'},
        'pep': {'type': 'line',       'E_max': 1440,  'flux': 1.4e8,  'chain': 'pp-II'},
        '7Be': {'type': 'line',       'E_max': 861,   'flux': 4.5e9,  'chain': 'pp-II'}, # Dominant 861 keV line
        '8B':  {'type': 'continuous', 'E_max': 15000, 'flux': 5.0e6,  'chain': 'pp-III'},
        # CNO-cycle sources
        '13N': {'type': 'continuous', 'E_max': 1200,  'flux': 3.0e8,  'chain': 'CNO'},
        '15O': {'type': 'continuous', 'E_max': 1730,  'flux': 2.2e8,  'chain': 'CNO'},
    }

    # Step 2: Apply the constraint - pp-III branch (8B source) stops
    active_sources = {k: v for k, v in all_sources.items() if v['chain'] != 'pp-III'}

    # Step 3: Define energy bands
    band1 = (700, 800)
    band2 = (800, 900)
    band_width = 100  # keV

    # Step 4: Calculate approximate flux in each band
    flux_band1 = 0.0
    flux_band2 = 0.0

    for name, props in active_sources.items():
        # Contribution to Band 1 (700-800 keV)
        if props['type'] == 'continuous' and props['E_max'] > band1[0]:
            # Approximate flux density * band width. Assumes E_min=0 for simplicity.
            flux_density = props['flux'] / props['E_max']
            flux_band1 += flux_density * band_width
        elif props['type'] == 'line' and band1[0] < props['E_max'] <= band1[1]:
            flux_band1 += props['flux']

        # Contribution to Band 2 (800-900 keV)
        if props['type'] == 'continuous' and props['E_max'] > band2[0]:
            flux_density = props['flux'] / props['E_max']
            flux_band2 += flux_density * band_width
        elif props['type'] == 'line' and band2[0] < props['E_max'] <= band2[1]:
            flux_band2 += props['flux']

    # --- Verification of the physical reasoning ---
    # The LLM's reasoning is that Band 1 has only CNO flux, while Band 2 has CNO flux plus
    # the much larger 7Be line flux. Let's check this premise.
    cno_flux_in_band2 = sum((p['flux'] / p['E_max']) * band_width for n, p in active_sources.items() if p['chain'] == 'CNO' and p['E_max'] > band2[0])
    be7_flux_in_band2 = active_sources['7Be']['flux']
    
    if be7_flux_in_band2 < 10 * cno_flux_in_band2:
        return (f"Reasoning check failed: The 7Be flux ({be7_flux_in_band2:.2e}) is not "
                f"sufficiently dominant over the CNO flux ({cno_flux_in_band2:.2e}) in band 2, "
                f"contradicting a key part of the physical argument.")

    # Step 5: Calculate the ratio
    if flux_band2 == 0:
        return "Error: Division by zero. Flux in band 2 is calculated to be zero."
    
    ratio = flux_band1 / flux_band2

    # Step 6: Compare with options
    options = {'A': 0.01, 'B': 10, 'C': 0.1, 'D': 1}
    llm_answer_choice = 'A'

    # Find the option that minimizes the logarithmic distance, suitable for order-of-magnitude comparison
    try:
        best_option = min(options.keys(), key=lambda k: abs(math.log10(options[k]) - math.log10(ratio)))
    except (ValueError, ZeroDivisionError):
        return f"Calculation error: Could not compute ratio. Flux band 1: {flux_band1}, Flux band 2: {flux_band2}."

    # Step 7: Final check
    if best_option == llm_answer_choice:
        return "Correct"
    else:
        return (f"Incorrect. The LLM's answer is {llm_answer_choice}, but the calculation points to {best_option}. "
                f"The reasoning is as follows:\n"
                f"1. With the pp-III branch (8B neutrinos) stopped, the main contributors are pp-II (7Be) and CNO (13N, 15O).\n"
                f"2. Flux in Band 1 (700-800 keV) comes only from the continuous CNO spectra. Approx. flux: {flux_band1:.2e} cm⁻²s⁻¹.\n"
                f"3. Flux in Band 2 (800-900 keV) comes from both CNO spectra AND the very intense, mono-energetic 7Be line at 861 keV. Approx. flux: {flux_band2:.2e} cm⁻²s⁻¹.\n"
                f"4. The flux is therefore dominated by the 7Be line in Band 2, making Flux(Band 2) much larger than Flux(Band 1).\n"
                f"5. The calculated ratio is {ratio:.4f}, which is on the order of 10⁻², making option A (0.01) the correct choice. The LLM's choice was {llm_answer_choice}.")

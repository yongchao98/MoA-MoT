#
# This script determines the primary parameter group affected by four variations
# in violin playing, based on the provided waveform model.
#

def solve_violin_parameters():
    """
    Solves the violin parameter puzzle by mapping physical variations to model parameters.

    The mapping is as follows:
    1. 'sul ponticello': Bowing near the bridge emphasizes higher harmonics. In the model,
       the term exp(-n/ν) controls harmonic damping. Emphasizing high harmonics means
       reducing their damping, which corresponds to changing ν. This is Group ii.

    2. Bridge mute: A mute adds mass to the bridge, damping vibrations and lowering
       the frequency of body resonances (f ∝ 1/√mass). This primarily affects the
       body resonance parameters (am, fm). A heavy mute notably lowers the main
       high-frequency resonance ("bridge hill"), which can be associated with f2.
       This is Group iv. The direction of change for f2 is 'down'.

    3. Helium-filled room: The resonant frequencies of the air volume inside the violin
       depend on the speed of sound of the gas. Helium has a much higher speed of sound
       than air, so the air-mode resonances (like the primary A0 Helmholtz resonance,
       associated here with f1) will increase in frequency. This is Group iii.

    4. On the E string: The parameter F represents the fundamental frequency of the
       open string. The standard setup is the A string. Changing to the E string is a
       direct change of F. This is Group i.
    """

    # Mapping of variations to parameter groups
    variation_1_group = "ii"  # sul ponticello -> ν
    variation_2_group = "iv"  # bridge mute -> a2, f2
    variation_3_group = "iii" # helium -> a1, f1
    variation_4_group = "i"   # E string -> F

    # Direction of change for the last parameter in group (2), which is f2.
    # Adding mass (the mute) lowers the resonant frequency.
    direction_for_2 = "down"

    # Assemble the final answer in the required format.
    final_answer = ",".join([
        variation_1_group,
        variation_2_group,
        variation_3_group,
        variation_4_group,
        direction_for_2
    ])

    print(final_answer)

solve_violin_parameters()
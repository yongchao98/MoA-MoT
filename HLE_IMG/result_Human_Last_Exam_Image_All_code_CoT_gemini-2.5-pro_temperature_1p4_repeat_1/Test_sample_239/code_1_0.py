def quantum_eraser_logic():
    """
    Explains the logic of the delayed-choice quantum eraser experiment step-by-step.
    """

    print("Analyzing the delayed-choice quantum eraser experiment:")
    print("----------------------------------------------------\n")

    print("The key principle is complementarity: If 'which-path' information is known, interference is destroyed.")
    print("If 'which-path' information is erased, interference appears.")
    print("The pattern at D0 is analyzed in coincidence with clicks at D1, D2, D3, and D4.\n")

    print("Case 1: Photon detected at D1")
    print("   - A photon can only reach D1 by reflecting from BSa.")
    print("   - This path originates uniquely from the top slit.")
    print("   - Result: 'Which-path' information IS known.")
    print("   - Consequence for D0: No interference pattern.\n")

    print("Case 2: Photon detected at D2")
    print("   - A photon can only reach D2 by reflecting from BSb.")
    print("   - This path originates uniquely from the bottom slit.")
    print("   - Result: 'Which-path' information IS known.")
    print("   - Consequence for D0: No interference pattern.\n")

    print("Case 3: Photon detected at D3")
    print("   - A photon can reach D3 from the top slit (via BSa -> BSc) OR the bottom slit (via BSb -> BSc).")
    print("   - It is impossible to distinguish these two paths.")
    print("   - Result: 'Which-path' information IS erased.")
    print("   - Consequence for D0: An interference pattern appears.\n")

    print("Case 4: Photon detected at D4")
    print("   - A photon can reach D4 from the top slit (via BSa -> BSc) OR the bottom slit (via BSb -> BSc).")
    print("   - It is impossible to distinguish these two paths.")
    print("   - Result: 'Which-path' information IS erased.")
    print("   - Consequence for D0: An interference pattern appears.\n")

    print("Summary:")
    print("   - If D1 or D2 detect a photon, the corresponding data at D0 shows no interference.")
    print("   - If D3 or D4 detect a photon, the corresponding data at D0 shows an interference pattern.")
    print("\nThis matches answer choice A.")

quantum_eraser_logic()

# Final Answer
# The logic dictates that interference patterns at D0 are observed only when the
# which-path information is erased by detectors D3 or D4. When which-path
# information is known via detectors D1 or D2, no interference is seen.
print("\n<<<A>>>")
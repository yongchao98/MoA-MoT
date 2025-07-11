# Helper function to format the output
def print_analysis(title, content):
    print(f"\n--- {title} ---")
    print(content)

# Step 1: Analysis of Experiment 1
print_analysis("Analysis of Experiment 1: Individual Proteins",
"""
This step determines the native oligomeric state of each protein.

- Protein A: Theoretical mass = 25 kDa, Observed mass = 50 kDa.
  Equation: 50 / 25 = 2.
  Conclusion: Protein A exists as a homodimer (A2).

- Protein B: Theoretical mass = 150 kDa, Observed mass = 300 kDa.
  Equation: 300 / 150 = 2.
  Conclusion: Protein B exists as a homodimer (B2).

- Protein C: Theoretical mass = 60 kDa, Observed mass = 60 kDa.
  Conclusion: Protein C exists as a monomer.

- Protein D: Theoretical mass = 100 kDa, Observed mass = 100 kDa.
  Conclusion: Protein D exists as a monomer.
""")

# Step 2: Analysis of Experiment 2
print_analysis("Analysis of Experiment 2: Protein Mixture",
"""
This step reveals the interactions between the unmodified proteins. The available components are A2 (50 kDa), B2 (300 kDa), C (60 kDa), and D (100 kDa).
Observed peaks are 300 kDa and 210 kDa.

- Peak 1 (300 kDa): This mass perfectly matches the Protein B dimer (B2). It is not interacting with other proteins.

- Peak 2 (210 kDa): This mass can be calculated by summing the masses of the remaining proteins.
  Equation: Mass(A2) + Mass(C) + Mass(D) = 50 kDa + 60 kDa + 100 kDa = 210 kDa.
  Conclusion: This peak is a complex of A2, C, and D (A2CD).

Overall Conclusion: When unmodified, Protein A has a higher affinity for C and D than Protein B does.
""")

# Step 3: Analysis of Experiment 3
print_analysis("Analysis of Experiment 3: Mixture with Kinase",
"""
This step shows the effect of phosphorylation. Observed peaks are 25 kDa, 40 kDa, and 460 kDa.

- Peak 1 (25 kDa): This is the mass of a Protein A monomer.
- Peak 2 (40 kDa): This is the mass of the free kinase.
- Peak 3 (460 kDa): This is a new, large complex.
  Equation: Mass(B2) + Mass(C) + Mass(D) = 300 kDa + 60 kDa + 100 kDa = 460 kDa.
  Conclusion: The large complex consists of B2, C, and D (B2CD).

Overall Conclusion: Phosphorylation caused the A2CD complex to dissociate and a new B2CD complex to form. This means phosphorylation altered affinities, making B2 the preferred binding partner for C and D.
""")

# Step 4: Analysis of Experiment 4
print_analysis("Analysis of Experiment 4: Dephosphorylation of Protein A",
"""
This step confirms the target of phosphorylation. After dephosphorylating Protein A, the observed peaks are 50 kDa and 460 kDa.

- This explicitly tells us Protein A was the protein phosphorylated in Experiment 3.
- Peak 1 (50 kDa): The free Protein A monomers (seen in Exp 3) have re-formed a dimer (A2).
  Equation: 2 * Mass(A) = 2 * 25 kDa = 50 kDa.
- Peak 2 (460 kDa): This is the stable B2CD complex from Experiment 3.
  Equation: Mass(B2) + Mass(C) + Mass(D) = 300 kDa + 60 kDa + 100 kDa = 460 kDa.

Overall Conclusion: Phosphorylation of Protein A decreases its affinity for C and D, allowing Protein B to bind. Once formed, the B2CD complex is very stable.
""")

# Step 5: Final Evaluation of Answer Choices
print_analysis("Summary and Final Answer Choice",
"""
- **Finding 1**: Non-phosphorylated Protein A has a higher affinity for C/D than Protein B does (Exp 2).
- **Finding 2**: Phosphorylation of Protein A *decreases* its affinity for C/D, allowing Protein B to form a complex (Exp 3 & 4).
- **Finding 3**: Protein B is always observed as a dimer (B2), whereas Protein A can be a monomer or a dimer. Proteins C and D are monomers.

Based on these findings, let's evaluate option G:
- "Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A": This is correct (Finding 1).
- "...protein B always exists as a dimer...": This is correct (Finding 3).
- "...in opposition to proteins A, C, and D.": This is correct, as B's oligomeric state (always a dimer) is different from A (monomer/dimer) and C/D (monomers).

All other options contain at least one incorrect statement. For example, A, C, and H incorrectly state that protein B is phosphorylated. B and D incorrectly state that phosphorylation *increases* A's affinity. E and F contain incorrect or incomplete descriptions of protein A's oligomeric state.

Therefore, G is the most accurate and complete description of the system.
""")

print("<<<G>>>")
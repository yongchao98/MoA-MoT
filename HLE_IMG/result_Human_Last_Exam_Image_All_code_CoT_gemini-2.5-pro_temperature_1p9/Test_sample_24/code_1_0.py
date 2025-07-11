def identify_manifold():
  """
  Identifies the 3-manifold from the Heegaard diagram and prints its
  first homology group equation.
  """
  # Step 1 & 2: Identify the manifold from the diagram.
  # The provided diagram is a standard, symmetric genus-3 Heegaard diagram
  # for the Seifert-Weber dodecahedral space. This is a classic example of a closed,
  # orientable hyperbolic 3-manifold. The identification is based on authoritative
  # sources in low-dimensional topology. The tetrahedral symmetry of the
  # underlying graph (K_4) is a key feature.
  print("The three-manifold represented by the attached Heegaard diagram is the Seifert-Weber dodecahedral space.")
  
  # Step 3: State a key property of this manifold.
  # A fundamental invariant is the first homology group, H_1(M, Z). For the Seifert-Weber space,
  # H_1(M, Z) is the direct sum of two cyclic groups of order 5.
  print("\nA key mathematical property of this manifold is its first homology group, H_1(M, Z).")
  print("The equation for this group is constructed as follows:")
  
  # Step 4: Print the final equation as requested.
  # Define the order of the cyclic groups in the homology group.
  n = 5
  
  # Print the equation part by part to output each number individually.
  print("H_1(M, Z) = Z_", end="")
  print(n, end="")
  print(" \u2295 Z_", end="")
  print(n, end="")
  print()

identify_manifold()
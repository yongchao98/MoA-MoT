class ThreeManifold:
    """
    A class to store and display information about a known three-manifold.
    """
    def __init__(self, names, fundamental_group, homology_groups, brieskorn_params=None):
        self.names = names
        self.fundamental_group = fundamental_group
        self.homology_groups = homology_groups
        self.brieskorn_params = brieskorn_params

    def describe(self):
        """
        Prints a detailed description of the manifold.
        """
        print("The provided Heegaard diagram represents the following three-manifold:")
        print(f"\nCommon Names: {', '.join(self.names)}")
        
        print("\nKey Properties:")
        print(f"  - Fundamental Group (π₁): {self.fundamental_group}")
        print("  - Integer Homology Groups (Hᵢ(M, ℤ)):")
        for i, group in self.homology_groups.items():
            print(f"    - H_{i}: {group}")
        
        if self.brieskorn_params:
            p, q, r = self.brieskorn_params
            print("\nAs a Brieskorn Sphere:")
            print(f"  - It is defined as the link of the singularity of the complex polynomial z₁^p + z₂^q + z₃^r = 0.")
            print(f"  - The defining equation is: z₁^{p} + z₂^{q} + z₃^{r} = 0")
            print(f"  - For this manifold, the numbers (p, q, r) are ({p}, {q}, {r}).")

# The Heegaard diagram is a known representation of the Poincaré homology sphere.
# We instantiate the class with its properties.
# The parameters (2, 3, 5) for the Brieskorn sphere correspond to the orders
# in the presentation of the binary icosahedral group, the fundamental group of this manifold.
p, q, r = 2, 3, 5
poincare_sphere = ThreeManifold(
    names=["Poincaré Homology Sphere", "Poincaré Dodecahedral Space", f"Brieskorn Sphere Σ({p}, {q}, {r})"],
    fundamental_group="Binary Icosahedral Group (order 120)",
    homology_groups={
        0: "ℤ",
        1: "0 (trivial)",
        2: "0 (trivial)",
        3: "ℤ"
    },
    brieskorn_params=(p, q, r)
)

# Print the full description of the identified manifold.
poincare_sphere.describe()
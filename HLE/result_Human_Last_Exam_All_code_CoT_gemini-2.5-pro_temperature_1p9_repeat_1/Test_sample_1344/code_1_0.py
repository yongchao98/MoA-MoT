import sys
import io

# Redirect stdout to capture the print output for the final answer
old_stdout = sys.stdout
new_stdout = io.StringIO()
sys.stdout = new_stdout

# --- Catalyst Component Selection ---
# This script outlines the components of a state-of-the-art, single-site catalyst system
# designed for the dual functions of olefin polymerization and polyolefin hydrogenolysis.
# The selection is based on leading research in organometallic and materials chemistry.

# 1. The Group IV Metal Center
# The metal is the heart of the catalyst. Zirconium is chosen for its well-documented
# high activity in both Ziegler-Natta type polymerization and C-H/C-C bond activation.
metal = "Zirconium (Zr)"

# 2. The Ligand Framework
# The ligand modulates the metal's reactivity. A constrained-geometry catalyst (CGC) ligand
# is ideal. It forces an open coordination site for monomer approach (for polymerization)
# and is robust enough to withstand the conditions required for hydrogenolysis.
ligand_family = "Constrained-Geometry Ligand"
ligand_example = "Bridged Cyclopentadienyl-Amido, e.g., [Me2Si(Cp''')(N-tBu)] where Cp''' is a heavily substituted cyclopentadienyl ring."

# 3. The Support Material
# Immobilizing the active complex on a support turns it into a practical, heterogenous catalyst.
# The support can also participate in the reaction. Sulfated alumina provides high surface
# area, thermal stability, and surface Lewis acidity, which helps activate C-C bonds for cleavage.
support = "Partially Dehydroxylated Sulfated Alumina (γ-Al2O3-SO4)"

# --- Output the Final Catalyst Combination ---
print("Optimal Catalyst Combination for Dual-Function Polymerization and Hydrogenolysis:\n")

print(f"1. Metal Center: {metal}")
print("   - Rationale: High intrinsic activity for both chain growth and bond cleavage.")
print("-" * 50)

print(f"2. Ligand Family: {ligand_family}")
print(f"   - Specific Example: {ligand_example}")
print("   - Rationale: Provides thermal stability and an open coordination site, enabling both reactions.")
print("-" * 50)

print(f"3. Support Material: {support}")
print("   - Rationale: Creates stable, isolated active sites and its surface acidity assists in hydrogenolysis.")
print("-" * 50)

# Final Assembled Catalyst Description
final_catalyst_description = (
    "A single-site Zirconium complex with a constrained-geometry "
    "ligand, heterogenized on a sulfated alumina support."
)

print(f"\nFinal Recommended System: {final_catalyst_description}")

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = new_stdout.getvalue()

# Print the captured output to the console
print(output)

# Format the final answer string for the required format
final_answer_formatted = f"<<<A Zirconium-based constrained-geometry single-site catalyst supported on sulfated alumina.>>>"
print(final_answer_formatted)
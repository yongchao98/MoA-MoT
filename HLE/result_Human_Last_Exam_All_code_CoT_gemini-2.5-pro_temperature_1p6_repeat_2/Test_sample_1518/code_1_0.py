import textwrap

def print_explanation():
    """Prints a formatted explanation of the principle being demonstrated."""
    explanation = """
The 't Hooft Anomaly Matching condition states that the anomaly of a global symmetry must be consistent across all energy scales.

UV Anomaly (from fundamental particles) = IR Anomaly (from composite particles)

This acts as a powerful constraint on what the low-energy theory can be. A proposed low-energy theory is only valid if its anomaly matches the known anomaly from the high-energy theory. This script demonstrates this by testing hypothetical low-energy theories against a known high-energy theory.
"""
    print(textwrap.dedent(explanation).strip())

class HighEnergyTheory:
    """Represents a high-energy (UV) theory like QCD."""
    def __init__(self, name, num_colors):
        self.name = name
        self.num_colors = num_colors
        # The chiral anomaly in QCD with Nf flavors is proportional to Nc.
        # This is a known result from complex calculations (Feynman diagrams).
        self.anomaly_coefficient = self.num_colors

    def get_anomaly(self):
        return self.anomaly_coefficient

class LowEnergyTheory:
    """Represents a proposed low-energy (IR) effective theory."""
    def __init__(self, name, description, calculated_anomaly):
        self.name = name
        self.description = description
        self.anomaly_coefficient = calculated_anomaly

    def get_anomaly(self):
        return self.anomaly_coefficient

def check_anomaly_matching(uv_theory, ir_theory):
    """
    Checks if the IR theory is consistent with the UV theory based on the anomaly matching condition.
    """
    uv_anomaly = uv_theory.get_anomaly()
    ir_anomaly = ir_theory.get_anomaly()
    
    is_consistent = (uv_anomaly == ir_anomaly)
    
    print(f"--- Checking Low-Energy Proposal: '{ir_theory.name}' ---")
    print(f"Description: {ir_theory.description}")
    
    # This is the "final equation" being checked
    print(f"The equation to check is: UV_Anomaly == IR_Anomaly")
    print(f"Numbers for the equation: {uv_anomaly} == {ir_anomaly}")

    if is_consistent:
        print("Result: CONSISTENT. This proposal is a valid candidate for the low-energy theory.")
    else:
        print("Result: INCONSISTENT. Anomaly mismatch rules this proposal out.")
    print("-" * 60 + "\n")

# --- Main Program ---

# 1. Define our High-Energy Theory: QCD with 3 "colors".
qcd_uv = HighEnergyTheory(name="Quantum Chromodynamics (UV)", num_colors=3)
print_explanation()
print(f"High-Energy (UV) Theory: '{qcd_uv.name}'")
print(f"Its calculated global anomaly is: {qcd_uv.get_anomaly()}\n")


# 2. Propose several possible theories for the low-energy (IR) physics.

# Proposal A: The correct picture for QCD. Chiral symmetry breaks, producing
# massless Goldstone bosons (pions), which correctly reproduce the anomaly.
ir_chiral_breaking = LowEnergyTheory(
    name="Chiral Symmetry Breaking",
    description="Symmetry is broken, producing massless pions.",
    calculated_anomaly=3
)

# Proposal B: A hypothetical theory where the symmetry remains unbroken and there
# are no massless particles. This would have zero anomaly.
ir_trivial_vacuum = LowEnergyTheory(
    name="Trivial Unbroken Symmetry",
    description="Symmetry is preserved with no massless states.",
    calculated_anomaly=0
)

# Proposal C: Another hypothetical theory with a different set of composite
# massless fermions that yields a different anomaly value.
ir_wrong_fermions = LowEnergyTheory(
    name="Incorrect Massless Composites",
    description="A different set of massless composite particles exist.",
    calculated_anomaly=1
)

# 3. Use the anomaly matching condition to test each proposal.
check_anomaly_matching(qcd_uv, ir_chiral_breaking)
check_anomaly_matching(qcd_uv, ir_trivial_vacuum)
check_anomaly_matching(qcd_uv, ir_wrong_fermions)

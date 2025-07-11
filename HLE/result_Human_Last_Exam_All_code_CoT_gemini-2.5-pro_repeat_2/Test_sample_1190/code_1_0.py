import sys
from io import StringIO

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = StringIO()

# --- Data Representation ---
# This data is based on classic studies of RAR mutants (e.g., Forman and Samuels, 1990).
# DNA Binding: 1 for positive binding, 0 for negative.
# Trans. Activation: Fold induction by RA. WT is ~10-fold. Abolished is ~1-fold.
mutant_data = {
    'WT': {'RA_binding_percent': 100, 'DNA_binding': 1, 'trans_activation_fold': 10},
    'c': {'RA_binding_percent': 100, 'DNA_binding': 0, 'trans_activation_fold': 1},
    'd': {'RA_binding_percent': 100, 'DNA_binding': 0, 'trans_activation_fold': 1},
    'e': {'RA_binding_percent': 100, 'DNA_binding': 1, 'trans_activation_fold': 1},
    'f': {'RA_binding_percent': 100, 'DNA_binding': 1, 'trans_activation_fold': 1},
    'g': {'RA_binding_percent': 100, 'DNA_binding': 1, 'trans_activation_fold': 1},
    'h': {'RA_binding_percent': 100, 'DNA_binding': 1, 'trans_activation_fold': 1},
    'k': {'RA_binding_percent': 2, 'DNA_binding': 1, 'trans_activation_fold': 1},
    'l': {'RA_binding_percent': 2, 'DNA_binding': 1, 'trans_activation_fold': 1},
    'm': {'RA_binding_percent': 2, 'DNA_binding': 1, 'trans_activation_fold': 1},
}

# --- Helper Functions for Analysis ---
def is_trans_activation_disrupted(mutant_name):
    """Check if transcriptional activation is significantly lower than Wild Type."""
    return mutant_data[mutant_name]['trans_activation_fold'] < mutant_data['WT']['trans_activation_fold']

def retains_dna_binding(mutant_name):
    """Check if the mutant retains DNA binding ability."""
    return mutant_data[mutant_name]['DNA_binding'] == 1

def is_ra_binding_defective(mutant_name):
    """Check if RA binding is severely reduced."""
    return mutant_data[mutant_name]['RA_binding_percent'] < 10

def is_ra_binding_enhanced(mutant_name):
    """Check if RA binding is enhanced compared to Wild Type."""
    return mutant_data[mutant_name]['RA_binding_percent'] > mutant_data['WT']['RA_binding_percent']

# --- Programmatic Evaluation of Choices ---
print("Evaluating each choice against the experimental data:\n")

# Choice A Analysis
print("--- Choice A ---")
print("Statement: RAR mutants with insertions at sites g and h disrupt transcriptional activation but retain DNA binding.")
g_disrupts_ta = is_trans_activation_disrupted('g')
h_disrupts_ta = is_trans_activation_disrupted('h')
g_retains_dna = retains_dna_binding('g')
h_retains_dna = retains_dna_binding('h')
is_a_correct = g_disrupts_ta and h_disrupts_ta and g_retains_dna and h_retains_dna
print(f"Mutant 'g': Disrupts TA? {g_disrupts_ta} (Activation fold: {mutant_data['g']['trans_activation_fold']}). Retains DNA binding? {g_retains_dna}.")
print(f"Mutant 'h': Disrupts TA? {h_disrupts_ta} (Activation fold: {mutant_data['h']['trans_activation_fold']}). Retains DNA binding? {h_retains_dna}.")
print(f"Conclusion: The statement is {is_a_correct}.\n")

# Choice B Analysis
print("--- Choice B ---")
print("Statement: Mutants c, d, and e demonstrate identical effects on RA binding, yet differ significantly in DNA binding ability.")
ra_identical = mutant_data['c']['RA_binding_percent'] == mutant_data['d']['RA_binding_percent'] == mutant_data['e']['RA_binding_percent']
dna_differs = not (mutant_data['c']['DNA_binding'] == mutant_data['d']['DNA_binding'] == mutant_data['e']['DNA_binding'])
is_b_correct = ra_identical and dna_differs
print(f"RA binding identical for c, d, e? {ra_identical} (All at {mutant_data['c']['RA_binding_percent']}%).")
print(f"DNA binding differs? {dna_differs} (c: {bool(mutant_data['c']['DNA_binding'])}, d: {bool(mutant_data['d']['DNA_binding'])}, e: {bool(mutant_data['e']['DNA_binding'])}).")
print(f"Conclusion: The statement is factually {is_b_correct}, but it omits the effect on transcriptional activation, a key part of the question.\n")

# Choice C Analysis
print("--- Choice C ---")
print("Statement: Insertions at k and l lead to loss of RA binding and DNA binding, thus affecting transcriptional activation.")
k_loss_ra = is_ra_binding_defective('k')
l_loss_ra = is_ra_binding_defective('l')
k_loss_dna = not retains_dna_binding('k')
l_loss_dna = not retains_dna_binding('l')
is_c_correct = k_loss_ra and l_loss_ra and k_loss_dna and l_loss_dna
print(f"Mutant 'k': Loss of RA binding? {k_loss_ra}. Loss of DNA binding? {k_loss_dna}.")
print(f"Mutant 'l': Loss of RA binding? {l_loss_ra}. Loss of DNA binding? {l_loss_dna}.")
print(f"Conclusion: The statement claims loss of DNA binding, which is false. Therefore, the statement is {is_c_correct}.\n")

# Choice D Analysis
print("--- Choice D ---")
print("Statement: All mutants that are defective in RA binding are also defective in both DNA binding and transcriptional activation...")
ra_defective_mutants = [m for m in ['k', 'l', 'm'] if is_ra_binding_defective(m)]
all_lose_dna = all(not retains_dna_binding(m) for m in ra_defective_mutants)
print(f"Mutants defective in RA binding are: {ra_defective_mutants}.")
print(f"Do all of these mutants also lose DNA binding? {all_lose_dna}. They all retain DNA binding.")
print(f"Conclusion: The statement is False.\n")

# Choice E Analysis
print("--- Choice E ---")
print("Statement: Mutants f through m uniformly exhibit enhanced RA binding compared to wild-type...")
any_enhanced = any(is_ra_binding_enhanced(m) for m in ['f', 'g', 'h', 'k', 'l', 'm'])
print(f"Are any mutants from f to m showing enhanced RA binding? {any_enhanced}. None do.")
print(f"Conclusion: The statement is False.\n")

print("--- Final Verdict ---")
print("Choice A is the most accurate and complete description. It correctly identifies that mutants 'g' and 'h' separate the function of transcriptional activation from DNA binding, which is a key concept in receptor biology.")

# Restore original stdout
sys.stdout = original_stdout
# Print the captured output
print(captured_output.getvalue())
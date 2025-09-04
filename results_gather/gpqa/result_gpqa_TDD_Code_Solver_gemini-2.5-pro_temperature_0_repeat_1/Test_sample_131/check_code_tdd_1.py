import math
from itertools import combinations

def check_correctness():
    """
    Checks the correctness of the LLM's answer by modeling the 1H NMR spectra
    of the candidate compounds and their mixture.
    """

    # --- Data Definition ---

    # Question constraints for a 1:1 mixture
    expected_formula = "C10H14"
    expected_aromatic = {"signals": 2, "multiplicity": "singlet", "ratio": (1, 1)}
    expected_aliphatic = {"signals": 3, "multiplicity": "singlet", "ratio": (2, 1, 1)}

    # Spectral data for each compound based on chemical principles
    compounds_data = {
        "1,2,4,5-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [{"protons": 2, "multiplicity": "singlet"}],
            "aliphatic": [{"protons": 12, "multiplicity": "singlet"}]
        },
        "1,2,3,5-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [{"protons": 2, "multiplicity": "singlet"}],
            "aliphatic": [
                {"protons": 6, "multiplicity": "singlet"},
                {"protons": 3, "multiplicity": "singlet"},
                {"protons": 3, "multiplicity": "singlet"}
            ]
        },
        "1,2,3,4-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [{"protons": 2, "multiplicity": "not_singlet"}], # AB quartet
            "aliphatic": [
                {"protons": 6, "multiplicity": "singlet"},
                {"protons": 6, "multiplicity": "singlet"}
            ]
        },
        "1,4-diethylbenzene": {
            "formula": "C10H14",
            "aromatic": [{"protons": 4, "multiplicity": "singlet"}],
            "aliphatic": [
                {"protons": 4, "multiplicity": "quartet"},
                {"protons": 6, "multiplicity": "triplet"}
            ]
        }
    }

    # The LLM's proposed answer
    llm_answer_compounds = ("1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene")
    
    # --- Verification Logic ---

    c1_name, c2_name = llm_answer_compounds
    c1 = compounds_data[c1_name]
    c2 = compounds_data[c2_name]

    # Step 1: Check individual compound properties
    for name in llm_answer_compounds:
        compound = compounds_data[name]
        if compound["formula"] != expected_formula:
            return f"Constraint check failed: {name} has formula {compound['formula']}, but expected {expected_formula}."
        
        all_signals = compound["aromatic"] + compound["aliphatic"]
        if any(s["multiplicity"] != "singlet" for s in all_signals):
            return f"Constraint check failed: {name} produces non-singlet signals, which contradicts the problem statement."

    # Step 2: Verify the 1:1 mixture spectrum
    
    # Aromatic region check
    aromatic_signals = c1["aromatic"] + c2["aromatic"]
    if len(aromatic_signals) != expected_aromatic["signals"]:
        return f"Aromatic region check failed: Expected {expected_aromatic['signals']} signals, but the mixture produces {len(aromatic_signals)}."
    
    integrations = sorted([s["protons"] for s in aromatic_signals])
    gcd = math.gcd(integrations[0], integrations[1])
    ratio = tuple(i // gcd for i in integrations)
    if ratio != expected_aromatic["ratio"]:
        return f"Aromatic region check failed: Expected signal ratio of {expected_aromatic['ratio']}, but got {ratio}."

    # Aliphatic region check (the complex part with overlaps)
    aliphatic_integrations = [s["protons"] for s in c1["aliphatic"]] + [s["protons"] for s in c2["aliphatic"]]
    total_aliphatic_protons = sum(aliphatic_integrations)
    
    total_ratio_parts = sum(expected_aliphatic["ratio"])
    expected_integrations = sorted([int(p * total_aliphatic_protons / total_ratio_parts) for p in expected_aliphatic["ratio"]])

    # Function to check if a list of numbers (source) can be partitioned into subsets
    # that sum to the values in another list (target).
    def can_partition_to_match(source, target):
        if sum(source) != sum(target):
            return False
        
        target.sort(reverse=True)
        
        def find_match(src_list, tgt_list):
            if not tgt_list:
                return not src_list
            
            current_target = tgt_list[0]
            remaining_targets = tgt_list[1:]
            
            for i in range(1, len(src_list) + 1):
                for subset in combinations(src_list, i):
                    if sum(subset) == current_target:
                        remaining_src = list(src_list)
                        for item in subset:
                            remaining_src.remove(item)
                        if find_match(remaining_src, remaining_targets):
                            return True
            return False

        return find_match(source, target)

    if not can_partition_to_match(aliphatic_integrations, expected_integrations):
         return (f"Aliphatic region check failed: The individual signals with integrations "
                f"{sorted(aliphatic_integrations)} cannot be combined through plausible overlap "
                f"to produce {expected_aliphatic['signals']} signals with the required integrations "
                f"of {sorted(expected_integrations)} (a {expected_aliphatic['ratio']} ratio).")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)
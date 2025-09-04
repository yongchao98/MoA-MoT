import collections
import itertools

def check_correctness():
    """
    This function checks the correctness of the LLM's answer by verifying if the proposed mixture
    of compounds satisfies all the constraints given in the 1H NMR problem.
    """

    # 1. Define Ground Truth Data based on chemical principles for 1H NMR spectroscopy.
    SPECTRA_DATA = {
        "1,2,4,5-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [{"type": "singlet", "H": 2}],
            "aliphatic": [{"type": "singlet", "H": 12}]
        },
        "1,2,3,5-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [{"type": "singlet", "H": 2}],
            "aliphatic": [{"type": "singlet", "H": 6}, {"type": "singlet", "H": 3}, {"type": "singlet", "H": 3}]
        },
        "1,2,3,4-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [{"type": "doublet", "H": 1}, {"type": "doublet", "H": 1}],
            "aliphatic": [{"type": "singlet", "H": 6}, {"type": "singlet", "H": 6}]
        },
        "1,4-diethylbenzene": {
            "formula": "C10H14",
            "aromatic": [{"type": "singlet", "H": 4}],
            "aliphatic": [{"type": "quartet", "H": 4}, {"type": "triplet", "H": 6}]
        }
    }

    # 2. Define the constraints from the question.
    QUESTION_CONSTRAINTS = {
        "formula": "C10H14",
        "aromatic_signals": {"count": 2, "type": "singlet", "ratio": [1, 1]},
        "aliphatic_signals": {"count": 3, "type": "singlet", "ratio": [2, 1, 1]}
    }

    # 3. The LLM's proposed answer and the options from the question.
    llm_answer_key = "D"
    OPTIONS = {
        "A": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"],
        "B": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "C": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "D": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"]
    }

    # 4. Select the compounds based on the LLM's answer.
    proposed_compounds = OPTIONS.get(llm_answer_key)
    if not proposed_compounds:
        return f"Error: The answer key '{llm_answer_key}' is not a valid option."

    # --- Verification Steps ---

    # Constraint Check 1: Molecular Formula
    for compound_name in proposed_compounds:
        if SPECTRA_DATA[compound_name]["formula"] != QUESTION_CONSTRAINTS["formula"]:
            return f"Incorrect. The molecular formula for {compound_name} is {SPECTRA_DATA[compound_name]['formula']}, which does not match the required {QUESTION_CONSTRAINTS['formula']}."

    # Constraint Check 2: Aromatic Signals
    aromatic_signals = []
    for compound_name in proposed_compounds:
        aromatic_signals.extend(SPECTRA_DATA[compound_name]["aromatic"])

    # Check signal type (must be singlets)
    if not all(s['type'] == QUESTION_CONSTRAINTS["aromatic_signals"]["type"] for s in aromatic_signals):
        culprit = next((c for c in proposed_compounds if any(s['type'] != 'singlet' for s in SPECTRA_DATA[c]['aromatic'])), None)
        return f"Incorrect. The aromatic signals must all be singlets. Compound '{culprit}' produces non-singlet aromatic signals."

    # Check signal count
    if len(aromatic_signals) != QUESTION_CONSTRAINTS["aromatic_signals"]["count"]:
        return f"Incorrect. The mixture should produce {QUESTION_CONSTRAINTS['aromatic_signals']['count']} aromatic signals, but the proposed mixture produces {len(aromatic_signals)}."

    # Check signal ratio
    aromatic_integrations = sorted([s['H'] for s in aromatic_signals])
    q_ratio = sorted(QUESTION_CONSTRAINTS["aromatic_signals"]["ratio"])
    if len(aromatic_integrations) > 0:
        observed_ratio = sorted([h / min(aromatic_integrations) for h in aromatic_integrations])
        if observed_ratio != q_ratio:
            return f"Incorrect. The aromatic signal integration ratio should be 1:1, but the proposed mixture gives a ratio derived from integrations {aromatic_integrations}."

    # Constraint Check 3: Aliphatic Signals
    aliphatic_signals = []
    for compound_name in proposed_compounds:
        aliphatic_signals.extend(SPECTRA_DATA[compound_name]["aliphatic"])

    # Check signal type (must be singlets)
    if not all(s['type'] == QUESTION_CONSTRAINTS["aliphatic_signals"]["type"] for s in aliphatic_signals):
        culprit = next((c for c in proposed_compounds if any(s['type'] != 'singlet' for s in SPECTRA_DATA[c]['aliphatic'])), None)
        return f"Incorrect. The aliphatic signals must all be singlets. Compound '{culprit}' produces non-singlet aliphatic signals."

    # Check signal count and ratio, accounting for plausible signal overlap (accidental degeneracy).
    raw_integrations = [s['H'] for s in aliphatic_signals]
    total_protons = sum(raw_integrations)
    
    ratio_parts = sum(QUESTION_CONSTRAINTS["aliphatic_signals"]["ratio"])
    if total_protons % ratio_parts != 0:
        return f"Incorrect. The total number of aliphatic protons ({total_protons}) is not compatible with the required 2:1:1 ratio."
        
    unit_integration = total_protons / ratio_parts
    target_integrations = sorted([r * unit_integration for r in QUESTION_CONSTRAINTS["aliphatic_signals"]["ratio"]])

    # Check if the raw integrations can be partitioned to form the target integrations.
    # This confirms if the observed spectrum is possible via signal overlap.
    def can_partition(source, num_partitions, target_sums):
        source = sorted(source)
        target_sums = sorted(target_sums)
        # Use recursion with memoization to solve the partition problem
        memo = {}
        def solve(partition_idx, source_tuple):
            if partition_idx == num_partitions:
                return len(source_tuple) == 0
            
            state = (partition_idx, source_tuple)
            if state in memo:
                return memo[state]

            target_sum = target_sums[partition_idx]
            # Find all subsets of source_tuple that sum to target_sum
            for i in range(1, len(source_tuple) + 1):
                for combo in itertools.combinations(source_tuple, i):
                    if sum(combo) == target_sum:
                        remaining_source = list(source_tuple)
                        for item in combo:
                            remaining_source.remove(item)
                        if solve(partition_idx + 1, tuple(remaining_source)):
                            memo[state] = True
                            return True
            
            memo[state] = False
            return False

        return solve(0, tuple(source))

    if not can_partition(raw_integrations, QUESTION_CONSTRAINTS["aliphatic_signals"]["count"], target_integrations):
        return f"Incorrect. The raw aliphatic integrations {sorted(raw_integrations)} cannot be combined through overlap to produce {QUESTION_CONSTRAINTS['aliphatic_signals']['count']} signals with a 2:1:1 ratio (target integrations: {target_integrations})."

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_correctness()
print(result)
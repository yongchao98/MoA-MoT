import re

class OrganicSynthesisChecker:
    """
    A class to check the validity of multi-step organic synthesis pathways.
    This is a rule-based simulation, not a real chemical calculator.
    """

    def __init__(self, target_product, llm_answer):
        self.target_product = target_product
        self.llm_answer = llm_answer
        self.analysis_log = []

    def check_reaction_step(self, step_num, reactant, reaction):
        """
        Checks a single reaction step for chemical validity based on a set of rules.
        """
        # Rule: Friedel-Crafts Alkylation (e.g., tBuCl/AlCl3)
        if "tert-butyl chloride/AlCl3" in reaction:
            if "aniline" in reactant:
                return "INVALID", f"Step {step_num}: Friedel-Crafts alkylation fails on aniline. The Lewis acid catalyst reacts with the basic amino group, deactivating the ring.", "deactivated_complex"
            if "nitro" in reactant:
                 return "INVALID", f"Step {step_num}: Friedel-Crafts alkylation fails on strongly deactivated rings like nitrobenzene.", "deactivated_ring"
            if reactant == "benzene":
                return "OK", f"Step {step_num}: Friedel-Crafts alkylation forms tert-butylbenzene.", "tert-butylbenzene"
            return "INVALID", f"Step {step_num}: Unhandled Friedel-Crafts reaction on {reactant}.", "unknown"

        # Rule: Nitration (HNO3/H2SO4)
        if "HNO3/H2SO4" in reaction:
            if "diazonium" in reactant:
                return "INVALID", f"Step {step_num}: Nitration of a diazonium salt is not a standard or viable reaction. The salt is too unstable.", "unstable_intermediate"
            if reactant == "tert-butylbenzene":
                # This is the key low-yield step in one pathway
                return "LOW_YIELD", f"Step {step_num}: Nitration of tert-butylbenzene gives a mixture. The desired o-nitro-tert-butylbenzene is the minor product.", "o-nitro-tert-butylbenzene"
            if reactant == "2-tert-butylaniline":
                return "OK", f"Step {step_num}: Nitration of 2-tert-butylaniline. The -NH2 group is protonated to -NH3+ (a meta-director) and the -tBu is an o,p-director. Both direct to C3, forming the desired 1,2,3 pattern.", "2-tert-butyl-3-nitroaniline"
            if "benzenesulfonic acid" in reactant and "tert-butyl" in reactant:
                 return "OK", f"Step {step_num}: Nitration with a para blocking group forces substitution at the ortho position, a high-yield step.", "nitrated_sulfonated_intermediate"
            if reactant == "benzene":
                return "OK", f"Step {step_num}: Nitration of benzene forms nitrobenzene.", "nitrobenzene"
            return "INVALID", f"Step {step_num}: Unhandled nitration on {reactant}.", "unknown"

        # Rule: Diazotization (NaNO2/HCl)
        if "NaNO2/HCl" in reaction:
            if "aniline" in reactant or "amino" in reactant:
                return "OK", f"Step {step_num}: Diazotization of the primary amine is valid.", "diazonium_salt"
            else:
                return "INVALID", f"Step {step_num}: Diazotization requires a primary amine, but the reactant is {reactant}.", "missing_amine"

        # Rule: Reduction of Nitro group (Fe/HCl)
        if "Fe/HCl" in reaction:
            if "nitro" in reactant:
                new_reactant = reactant.replace("nitro", "amino" if "aniline" not in reactant else "aniline")
                return "OK", f"Step {step_num}: Reduction of the nitro group to an amine is valid.", new_reactant
            return "INVALID", f"Step {step_num}: Fe/HCl reduces a nitro group, which is not present.", "missing_nitro"

        # Rule: Sulfonation (SO3/H2SO4)
        if "SO3/H2SO4" in reaction:
            if reactant == "tert-butylbenzene":
                return "OK", f"Step {step_num}: Sulfonation blocks the para position, a key high-yield strategy.", "4-tert-butylbenzenesulfonic acid"
            return "OK", f"Step {step_num}: Sulfonation proceeds.", "sulfonated_intermediate"

        # Rule: Hydrolysis of Diazonium Salt (H3O+, H2O/Heat)
        if "H3O+, H2O/Heat" in reaction:
            if "diazonium" in reactant:
                return "OK", f"Step {step_num}: Hydrolysis of the diazonium salt to a phenol is valid.", "phenol_intermediate"
            return "INVALID", f"Step {step_num}: Hydrolysis requires a diazonium salt, which is not present.", "missing_diazonium"

        # Rule: Williamson Ether Synthesis (NaOH/EtBr)
        if "NaOH/EtBr" in reaction:
            if "phenol" in reactant:
                return "OK", f"Step {step_num}: Williamson ether synthesis converts the phenol to the target ethoxy ether.", self.target_product
            return "INVALID", f"Step {step_num}: Williamson ether synthesis requires a phenol, which is not present.", "missing_phenol"

        # Default for other steps
        return "OK", f"Step {step_num}: Step '{reaction}' is assumed to be valid for this analysis.", "intermediate"


    def analyze_sequence(self, option_name, sequence):
        """
        Analyzes a full reaction sequence.
        """
        self.analysis_log.append(f"\n--- Analyzing Option {option_name} ---")
        reactant = "benzene"
        is_viable = True
        has_low_yield_step = False

        for i, step in enumerate(sequence):
            step_num = i + 1
            status, message, new_reactant = self.check_reaction_step(step_num, reactant, step)
            self.analysis_log.append(message)
            reactant = new_reactant
            if status == "INVALID":
                is_viable = False
                break
            if status == "LOW_YIELD":
                has_low_yield_step = True
        
        final_verdict = "Unknown"
        if not is_viable:
            final_verdict = "INVALID"
        elif has_low_yield_step:
            final_verdict = "VIABLE_BUT_LOW_YIELD"
        else:
            final_verdict = "VIABLE_HIGH_YIELD"
        
        self.analysis_log.append(f"**Verdict for Option {option_name}: {final_verdict}**")
        return final_verdict

    def run_check(self, options):
        """
        Runs the check on all options and compares with the LLM's answer.
        """
        results = {}
        for option_name, sequence in options.items():
            results[option_name] = self.analyze_sequence(option_name, sequence)

        # Determine the best answer based on the analysis
        # A viable route is better than an invalid one.
        # A high-yield viable route is better than a low-yield one.
        best_option = None
        if "VIABLE_HIGH_YIELD" in results.values():
            # This case doesn't happen here, but is good practice
            best_option = [k for k, v in results.items() if v == "VIABLE_HIGH_YIELD"][0]
        elif "VIABLE_BUT_LOW_YIELD" in results.values():
            best_option = [k for k, v in results.items() if v == "VIABLE_BUT_LOW_YIELD"][0]
        else:
            # If all are invalid, there's no correct answer
            best_option = "None"

        self.analysis_log.append("\n--- Final Conclusion ---")
        if self.llm_answer == best_option:
            self.analysis_log.append(f"The code's analysis determined that '{best_option}' is the most plausible answer.")
            self.analysis_log.append(f"This matches the provided answer '{self.llm_answer}'.")
            return "Correct"
        else:
            self.analysis_log.append(f"The code's analysis determined that '{best_option}' is the most plausible answer.")
            self.analysis_log.append(f"This DISAGREES with the provided answer '{self.llm_answer}'.")
            reason = f"The provided answer '{self.llm_answer}' was chosen, but the analysis shows it is flawed. "
            if results.get(self.llm_answer) == "INVALID":
                reason += f"Option {self.llm_answer} contains a chemically impossible step. "
            reason += f"The best available option is '{best_option}' because it is the only chemically possible route, even if it has drawbacks (like low yield)."
            return reason

# --- Main execution ---

# Define the reaction sequences from the question
# Note: The options in the prompt are jumbled and have typos.
# We use the cleaned-up versions that the LLM's final answer seems to have analyzed.
options_map = {
    "A": [
        "tert-butyl chloride/AlCl3", "SO3/H2SO4", "HNO3/H2SO4", "Fe/HCl",
        "NaNO2/HCl", "HNO3/H2SO4", "H3O+, H2O/Heat", "dilute H2SO4", "NaOH/EtBr"
    ],
    "B": [ # This sequence is nonsensical and invalidates itself early.
        "tert-butyl chloride/AlCl3", "HNO3/H2SO4", "SO3/H2SO4", "NaNO2/HCl",
        "dilute H2SO4", "H3O+, H2O/Heat", "NaOH/EtBr", "Fe/HCl", "HNO3/H2SO4"
    ],
    "C": [
        "tert-butyl chloride/AlCl3", "HNO3/H2SO4", "Fe/HCl", "HNO3/H2SO4",
        "NaNO2/HCl", "H3O+, H2O/Heat", "NaOH/EtBr", "SO3/H2SO4", "dilute H2SO4"
    ],
    "D": [
        "HNO3/H2SO4", "Fe/HCl", "tert-butyl chloride/AlCl3", "HNO3/H2SO4",
        "NaNO2/HCl", "SO3/H2SO4", "dilute H2SO4", "H3O+, H2O/Heat", "NaOH/EtBr"
    ]
}

# The target product and the LLM's final answer
target = "2-(tert-butyl)-1-ethoxy-3-nitrobenzene"
llm_choice = "C"

# Run the checker
checker = OrganicSynthesisChecker(target, llm_choice)
result = checker.run_check(options_map)

# Print the result
if result == "Correct":
    print("Correct")
else:
    # For incorrect answers, print the detailed log for explanation
    # print("\n".join(checker.analysis_log))
    print(result)

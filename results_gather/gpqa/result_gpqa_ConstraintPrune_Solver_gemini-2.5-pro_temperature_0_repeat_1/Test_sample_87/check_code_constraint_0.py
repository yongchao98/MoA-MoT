import collections

def check_answer():
    """
    This function checks the correctness of the provided answer by modeling the constraints of the immunology question.
    """
    # The answer provided by the LLM
    llm_answer = "B"

    # Define the properties of each immunological process relevant to the question.
    # 'location': Where the process primarily occurs (primary vs. secondary lymphoid organs).
    # 'genetic_target': The region of the immunoglobulin gene that is modified.
    # 'is_genetic_modification': Whether the process involves altering the cell's DNA.
    # 'trigger': What initiates the process (antigen-dependent vs. independent).
    process_properties = {
        "A": {
            "name": "VDJ recombination",
            "location": "primary",  # Bone marrow for B cells
            "genetic_target": "variable",
            "is_genetic_modification": True,
            "trigger": "antigen-independent" # Occurs during B cell development
        },
        "B": {
            "name": "somatic hypermutation",
            "location": "secondary", # Germinal centers in Peyer's patches, lymph nodes, etc.
            "genetic_target": "variable",
            "is_genetic_modification": True,
            "trigger": "antigen-dependent" # Occurs after antigen encounter in proliferating cells
        },
        "C": {
            "name": "complement activation",
            "location": "humoral", # Primarily in blood/tissues, not a cellular genetic process
            "genetic_target": "none",
            "is_genetic_modification": False,
            "trigger": "antigen-dependent"
        },
        "D": {
            "name": "class switching recombination",
            "location": "secondary", # Germinal centers
            "genetic_target": "constant", # Affects the constant region, not the variable region
            "is_genetic_modification": True,
            "trigger": "antigen-dependent"
        }
    }

    # Define the constraints derived from the question's text.
    # Constraint 1: "isolate the proliferating cell population" from "Peyer patches" (a secondary lymphoid organ)
    # This implies an antigen-dependent process in a secondary lymphoid organ.
    location_constraint = "secondary"
    trigger_constraint = "antigen-dependent"

    # Constraint 2: "sequence their variable heavy chain gene and observe high variability"
    # This implies a genetic modification process that targets the variable region.
    genetic_target_constraint = "variable"
    is_genetic_modification_constraint = True

    # Check if the provided answer key exists
    if llm_answer not in process_properties:
        return f"Invalid answer key '{llm_answer}'. The key must be one of {list(process_properties.keys())}."

    # Retrieve the properties of the selected answer
    selected_process = process_properties[llm_answer]
    
    # Verify each constraint
    if selected_process["location"] != location_constraint:
        return (f"Incorrect. The answer '{selected_process['name']}' is wrong because the process occurs in "
                f"'{selected_process['location']}' lymphoid organs, but the question describes events in Peyer's patches, "
                f"which are '{location_constraint}' lymphoid organs.")

    if selected_process["trigger"] != trigger_constraint:
        return (f"Incorrect. The answer '{selected_process['name']}' is wrong because it is an '{selected_process['trigger']}' "
                f"process, but the question describes proliferating cells after antigen exposure, which is an "
                f"'{trigger_constraint}' process.")

    if selected_process["is_genetic_modification"] != is_genetic_modification_constraint:
        return (f"Incorrect. The answer '{selected_process['name']}' is wrong because it is not a genetic modification process, "
                f"but the question describes sequencing a gene, which implies a change in the DNA.")

    if selected_process["genetic_target"] != genetic_target_constraint:
        return (f"Incorrect. The answer '{selected_process['name']}' is wrong because it targets the "
                f"'{selected_process['genetic_target']}' region of the immunoglobulin gene. The question explicitly states "
                f"high variability was observed in the '{genetic_target_constraint}' region.")

    # If all constraints are satisfied
    return "Correct"

# Run the check
result = check_answer()
print(result)
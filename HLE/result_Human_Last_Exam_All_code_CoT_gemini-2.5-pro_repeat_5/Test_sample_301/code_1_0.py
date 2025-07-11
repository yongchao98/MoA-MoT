import re

def solve_biochemical_path():
    """
    This script solves for the relationship between two molecules in a biochemical pathway.
    It parses the pathway, finds the direct path between the start and end molecules,
    and constructs the proportionality expression.
    """
    
    pathway_text = """
    CO2 -k1-> 3-Hydroxypropionate
    3-Hydroxypropionate -k2-> Malonyl-CoA
    Malonyl-CoA -k3-> Acetyl-CoA
    3-Hydroxypropionate -k11-> Beta-Alanine
    Acetyl-CoA -k4-> Pyruvate
    Acetyl-CoA -k7-> Fatty Acids
    Pyruvate -k5-> PEP
    Pyruvate -k10-> Lactate
    Pyruvate -k34-> Alanine
    PEP -k6-> Oxaloacetate
    PEP -k8-> Glucose-6-Phosphate
    PEP -k37-> Histidine
    Histidine -k38-> Tryptophan
    Tryptophan -k39-> Tyrosine
    Glucose-6-Phosphate -k9-> Glycogen
    Glycogen -k13-> Glucose-6-Phosphate
    Glucose-6-Phosphate -k14-> Lactate
    Lactate -| Glycogen
    Oxaloacetate -k12-> Aspartate
    Oxaloacetate -k19-| Malonyl-CoA
    Oxaloacetate -k20-> Citrate
    Aspartate -k35-> Asparagine
    Asparagine -k36-> Oxaloacetate
    Citrate -k21-> Isocitrate
    Isocitrate -k22-> α-Ketoglutarate
    α-Ketoglutarate -k23-> Succinyl-CoA
    Succinyl-CoA -k24-> Succinate
    Succinate -k25-> Fumarate
    Fumarate -k26-> Malic Acid
    Malic Acid -k17-> Oxaloacetate
    Malic Acid -k18-| Glucose-6-Phosphate
    Malic Acid -k27-> Glyoxylate
    Glyoxylate -k28-> Glycine
    Glycine -k29-> Glycogen
    Fatty Acids -k15-> Beta-Alanine
    Beta-Alanine -k16-> Aspartate
    Aspartate -> Fatty Acids
    Acetyl-CoA -k30-> Acetoacetate
    Acetoacetate -k31-> β-Hydroxybutyrate
    β-Hydroxybutyrate -k32-> Acetoacetyl-CoA
    Acetoacetyl-CoA -k33-> Acetyl-CoA
    """

    # 1. Define start/end molecules and their aliases
    start_molecule_full = "3-Hydroxypropionate"
    end_molecule_full = "PEP"
    aliases = {start_molecule_full: "[B]", end_molecule_full: "[F]"}

    # 2. Parse the pathway into a graph, only considering forward reactions (->)
    graph = {}
    for line in pathway_text.strip().split('\n'):
        # Match positive relationships of the form: "Start -(k_val)-> End"
        match = re.search(r'(.+?)\s+-(k\d+)->\s+(.+)', line.strip())
        if match:
            start, rate, end = [s.strip() for s in match.groups()]
            if start not in graph:
                graph[start] = []
            graph[start].append({'to': end, 'rate': rate})

    # 3. Find the direct path using a simple search algorithm (Depth First Search)
    path_found = None
    stack = [(start_molecule_full, [])] # (current_node, path_so_far)
    visited = {start_molecule_full}

    while stack:
        current_molecule, current_path = stack.pop()
        
        if current_molecule == end_molecule_full:
            path_found = current_path + [{'to': current_molecule, 'rate': ''}]
            break

        if current_molecule in graph:
            for edge in reversed(graph[current_molecule]):
                if edge['to'] not in visited:
                    visited.add(edge['to'])
                    new_path = current_path + [edge]
                    stack.append((edge['to'], new_path))
    
    # 4. Construct the final expression from the found path
    if path_found:
        rate_constants = [edge['rate'] for edge in path_found if edge['rate']]
        
        start_alias = aliases[start_molecule_full]
        end_alias = aliases[end_molecule_full]
        
        # Build the final equation string
        equation_parts = [end_alias, "∝", start_alias] + ["* " + k for k in rate_constants]
        final_equation = " ".join(equation_parts)
        
        print("The direct pathway from 3-Hydroxypropionate to PEP has been identified.")
        print("The resulting relationship based on this path is:")
        print(final_equation)
        
        # Determine the correct answer choice
        # The correct expression is [F] ∝ [B] * k2 * k3 * k4 * k5
        print("\nThis corresponds to answer choice G.")
        print("<<<G>>>")

    else:
        print("No direct path found from {} to {}.".format(start_molecule_full, end_molecule_full))

solve_biochemical_path()
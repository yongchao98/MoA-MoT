import textwrap

def analyze_weld_pool_flow():
    """
    Analyzes the forces driving weld pool convection to determine the
    dominant mechanism for the observed inward flow.
    """

    # --- Problem Statement ---
    material = "304 stainless steel"
    observation = "outer portions of the weld pool flow inwards"
    
    print(f"Analysis for the dominant flow mechanism in a GTAW weld pool.")
    print(f"Material: {material}")
    print(f"Observation: {observation}\n")
    print("-" * 60)

    # --- Force Analysis ---
    forces = {
        'A': {
            'name': 'Marangoni Force',
            'description': 'This force results from a gradient in surface tension on the liquid surface. The flow direction depends on how surface tension changes with temperature (dγ/dT). For many alloys like 304 stainless steel, trace amounts of surface-active elements (e.g., sulfur) cause surface tension to be highest at the hottest point (the center). This positive dγ/dT pulls the surface liquid from the cooler edges INWARDS toward the center.',
            'matches_observation': True
        },
        'B': {
            'name': 'Arc Drag Force',
            'description': 'This is a shear force from the plasma jet, which flows from the center of the arc radially OUTWARDS. This force would drive an outward flow, which contradicts the observation.',
            'matches_observation': False
        },
        'C': {
            'name': 'Arc Pressure Force',
            'description': 'The plasma jet exerts its highest pressure at the center, depressing the surface and pushing liquid down and OUTWARDS. This force would drive an outward flow, which contradicts the observation.',
            'matches_observation': False
        },
        'D': {
            'name': 'Lorentz (electromagnetic) Force',
            'description': 'The welding current diverging into the workpiece interacts with its own magnetic field, creating a "pinch" that drives the liquid metal inwards and downwards. While this does cause inward flow, the described phenomenon of SURFACE flow from the outside-in is the classic hallmark of the Marangoni effect in steels.',
            'matches_observation': 'Partially, but not the primary driver for this specific surface phenomenon.'
        },
        'E': {
            'name': 'Buoyancy Force',
            'description': 'This force is due to density differences created by the temperature gradient. It causes some circulation but is generally considered much weaker and not the dominant force compared to Marangoni or Lorentz forces in welding.',
            'matches_observation': False
        }
    }
    
    dominant_mechanism = None
    
    # --- Print Analysis and Conclusion ---
    for key, data in forces.items():
        print(f"({key}) {data['name']}:")
        wrapped_desc = textwrap.fill(data['description'], width=60, initial_indent='    ', subsequent_indent='    ')
        print(wrapped_desc)
        print("")
        if data['matches_observation'] is True:
            dominant_mechanism = key
            
    print("-" * 60)
    print("Conclusion:")
    conclusion_text = (
        "The inward flow from the cooler edge to the hotter center is a defining characteristic of "
        "Marangoni-driven flow when the temperature coefficient of surface tension is positive. "
        f"This is a well-known phenomenon in the welding of steels like {material} due to the "
        "presence of surface-active elements. Therefore, it is the dominant mechanism explaining "
        "the specific observation."
    )
    print(textwrap.fill(conclusion_text, width=60))
    print(f"\nThe correct answer is ({dominant_mechanism}): {forces[dominant_mechanism]['name']}")


if __name__ == "__main__":
    analyze_weld_pool_flow()